import pytest

from plasmidlab.core import (
    CloningError,
    Feature,
    MoleculeTopology,
    SequenceRecord,
    gibson_assembly,
    golden_gate_assembly,
    inverse_pcr_mutagenesis,
    restriction_clone,
)


def test_restriction_cloning_ligates_compatible_directional_ends_and_features() -> None:
    vector = SequenceRecord(
        id="vector",
        sequence="TTTTGAATTCGGGGGGATCCCCCC",
        topology=MoleculeTopology.CIRCULAR,
        features=(Feature(type="rep_origin", start=20, end=24, name="ori"),),
    )
    insert = SequenceRecord(
        id="insert",
        sequence="GAATTCATGAAAGGATCC",
        features=(Feature(type="CDS", start=6, end=12, name="payload"),),
    )

    result = restriction_clone(vector, insert, "EcoRI,BamHI")

    assert result.product.sequence == "GATCCCCCCTTTTGAATTCATGAAAG"
    assert result.product.topology is MoleculeTopology.CIRCULAR
    assert [feature.name for feature in result.product.features] == ["ori", "payload"]
    assert [(feature.start, feature.end) for feature in result.product.features] == [(5, 9), (19, 25)]
    assert result.product.history[-1].operation == "restriction_clone"
    assert result.product.history[-1].parameters["vector_enzymes"] == ("EcoRI", "BamHI")


def test_restriction_cloning_rejects_incompatible_ends() -> None:
    vector = SequenceRecord(
        id="vector",
        sequence="TTTTGAATTCGGGGGGATCCCCCC",
        topology=MoleculeTopology.CIRCULAR,
    )
    insert = SequenceRecord(id="bad_insert", sequence="GAATTCATGAAAAAGCTT")

    with pytest.raises(CloningError, match="incompatible ends"):
        restriction_clone(vector, insert, "EcoRI,BamHI", "EcoRI,HindIII")


def test_restriction_cloning_supports_reverse_insert_orientation() -> None:
    vector = SequenceRecord(
        id="single_cut_vector",
        sequence="AAAAGAATTCCCC",
        topology=MoleculeTopology.CIRCULAR,
    )
    insert = SequenceRecord(id="single_cut_insert", sequence="GAATTCAAAGGGGAATTC")

    result = restriction_clone(vector, insert, "EcoRI", insert_orientation="reverse")

    assert result.product.sequence == "AATTCCCCAAAAGCCCCTTTGAATT"
    assert result.product.history[-1].parameters["insert_orientation"] == "reverse"


def test_gibson_assembly_uses_20_to_40_bp_overlaps_exactly() -> None:
    overlap_a = "G" * 20
    overlap_b = "TACG" * 8
    fragments = (
        SequenceRecord(id="g1", sequence="AAAACCCC" + overlap_a),
        SequenceRecord(id="g2", sequence=overlap_a + "TTTTGGGG" + overlap_b),
        SequenceRecord(id="g3", sequence=overlap_b + "CCCCAAAA"),
    )

    result = gibson_assembly(fragments, min_overlap=20)

    assert result.product.sequence == "AAAACCCC" + overlap_a + "TTTTGGGG" + overlap_b + "CCCCAAAA"
    assert result.product.topology is MoleculeTopology.LINEAR
    assert [overlap.length for overlap in result.overlaps] == [20, 32]
    assert result.product.history[-1].operation == "gibson_assembly"


def test_gibson_circular_assembly_removes_terminal_overlap() -> None:
    overlap_a = "A" * 20
    overlap_b = "C" * 20
    overlap_c = "G" * 20
    fragments = (
        SequenceRecord(id="c1", sequence=overlap_c + "TTTT" + overlap_a),
        SequenceRecord(id="c2", sequence=overlap_a + "GGGG" + overlap_b),
        SequenceRecord(id="c3", sequence=overlap_b + "CCCC" + overlap_c),
    )

    result = gibson_assembly(fragments, min_overlap=20, circular=True)

    assert result.product.sequence == overlap_c + "TTTT" + overlap_a + "GGGG" + overlap_b + "CCCC"
    assert result.product.topology is MoleculeTopology.CIRCULAR
    assert [overlap.length for overlap in result.overlaps] == [20, 20, 20]


def test_golden_gate_bsa_i_like_type_iis_assembly() -> None:
    fragments = (
        _golden_gate_fragment("gg1", "AAAA", "CCCC", "TTTT"),
        _golden_gate_fragment("gg2", "TTTT", "GGGG", "AAAA"),
    )

    result = golden_gate_assembly(fragments, "BsaI", circular=True)

    assert result.product.sequence == "AAAACCCCTTTTGGGG"
    assert result.product.topology is MoleculeTopology.CIRCULAR
    assert result.warnings == ()
    assert [overlap.sequence for overlap in result.overlaps] == ["TTTT", "AAAA"]
    assert result.product.history[-1].operation == "golden_gate_assembly"


def test_golden_gate_uses_non_palindromic_sequence_derived_overhangs_exactly() -> None:
    fragments = (
        _golden_gate_fragment("gg1", "AAAA", "CCCC", "ATGC"),
        _golden_gate_fragment("gg2", "ATGC", "GGGG", "TTAA"),
        _golden_gate_fragment("gg3", "TTAA", "TTTT", "AAAA"),
    )

    result = golden_gate_assembly(fragments, "BsaI", circular=True)

    assert result.product.sequence == "AAAACCCCATGCGGGGTTAATTTT"
    assert [overlap.sequence for overlap in result.overlaps] == ["ATGC", "TTAA", "AAAA"]


def test_golden_gate_warns_on_duplicated_overhangs() -> None:
    fragments = (
        _golden_gate_fragment("gg1", "AAAA", "CCCC", "TTTT"),
        _golden_gate_fragment("gg2", "TTTT", "GGGG", "CCCC"),
        _golden_gate_fragment("gg3", "CCCC", "TTTT", "TTTT"),
        _golden_gate_fragment("gg4", "TTTT", "AAAA", "AAAA"),
    )

    result = golden_gate_assembly(fragments, "BsaI", circular=True)

    assert "duplicated left overhangs: TTTT" in result.warnings
    assert "duplicated right overhangs: TTTT" in result.warnings
    assert "ambiguous Golden Gate overhang reuse: TTTT" in result.warnings


def test_golden_gate_rejects_incompatible_overhang_order() -> None:
    fragments = (
        _golden_gate_fragment("gg1", "AAAA", "CCCC", "TTTT"),
        _golden_gate_fragment("gg2", "CCCC", "GGGG", "AAAA"),
    )

    with pytest.raises(CloningError, match="does not match"):
        golden_gate_assembly(fragments, "BsaI", circular=True)


def test_golden_gate_rejects_reversed_non_palindromic_fragment_orientation() -> None:
    first = _golden_gate_fragment("gg1", "AAAA", "CCCC", "ATGC")
    second = _golden_gate_fragment("gg2", "ATGC", "GGGG", "AAAA").reverse_complement()

    with pytest.raises(CloningError, match="does not match"):
        golden_gate_assembly((first, second), "BsaI", circular=True)


def test_golden_gate_rejects_ambiguous_overhangs_by_default() -> None:
    fragments = (
        _golden_gate_fragment("gg1", "AAAA", "CCCC", "ANNN"),
        _golden_gate_fragment("gg2", "ATGC", "GGGG", "AAAA"),
    )

    with pytest.raises(CloningError, match="does not match"):
        golden_gate_assembly(fragments, "BsaI", circular=True)


def test_inverse_pcr_mutagenesis_deletes_inserts_circularizes_and_tracks_provenance() -> None:
    template = SequenceRecord(
        id="template",
        sequence="AAAACCCCGGGGTTTT",
        topology=MoleculeTopology.CIRCULAR,
        features=(Feature(type="misc_feature", start=8, end=12, name="kept"),),
    )

    result = inverse_pcr_mutagenesis(template, 4, 8, insertion_sequence="GG")

    assert result.product.sequence == "AAAAGGGGGGTTTT"
    assert result.product.topology is MoleculeTopology.CIRCULAR
    assert result.product.features[0].start == 6
    assert result.product.features[0].end == 10
    assert result.product.history[-1].operation == "inverse_pcr_mutagenesis"
    assert result.product.history[-1].parameters["inserted_length"] == 2


def _golden_gate_fragment(
    record_id: str,
    left_overhang: str,
    payload: str,
    right_overhang: str,
) -> SequenceRecord:
    return SequenceRecord(
        id=record_id,
        sequence="GGTCTC" + "A" + left_overhang + payload + right_overhang + "A" + "GAGACC",
    )
