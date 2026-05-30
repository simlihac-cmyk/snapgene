from plasmidlab.core import (
    MoleculeTopology,
    SequenceRecord,
    gibson_assembly,
    golden_gate_assembly,
    restriction_clone,
    simulate_pcr,
)


def test_restriction_clone_and_gibson_events_include_full_provenance_fields() -> None:
    vector = SequenceRecord(
        id="vector",
        sequence="TTTTGAATTCGGGGGGATCCCCCC",
        topology=MoleculeTopology.CIRCULAR,
    )
    insert = SequenceRecord(id="insert", sequence="GAATTCATGAAAGGATCC")

    cloned = restriction_clone(vector, insert, "EcoRI,BamHI").product

    assert cloned.history[-1].operation == "restriction_clone"
    assert cloned.history[-1].input_record_ids == ("vector", "insert")
    assert cloned.history[-1].input_version_ids == (vector.version_id, insert.version_id)
    assert cloned.history[-1].output_record_id == cloned.id
    assert cloned.history[-1].output_version_id == cloned.version_id
    assert cloned.history[-1].affected_ranges == ((0, cloned.length),)
    assert "Restriction cloned" in cloned.history[-1].description

    overlap = "A" * 20
    fragments = (
        SequenceRecord(id="g1", sequence="CCCC" + overlap),
        SequenceRecord(id="g2", sequence=overlap + "GGGG"),
    )
    assembled = gibson_assembly(fragments, min_overlap=20).product

    assert assembled.history[-1].operation == "gibson_assembly"
    assert assembled.history[-1].input_record_ids == ("g1", "g2")
    assert assembled.history[-1].input_version_ids == tuple(fragment.version_id for fragment in fragments)
    assert assembled.history[-1].output_record_id == assembled.id
    assert assembled.history[-1].output_version_id == assembled.version_id
    assert assembled.history[-1].affected_ranges == ((0, assembled.length),)


def test_golden_gate_product_records_output_version_id() -> None:
    fragments = (
        _golden_gate_fragment("gg1", "AAAA", "CCCC", "TTTT"),
        _golden_gate_fragment("gg2", "TTTT", "GGGG", "AAAA"),
    )

    product = golden_gate_assembly(fragments, "BsaI").product
    event = product.history[-1]

    assert product.version_id.startswith("golden_gate_assembly:v")
    assert product.version_id != "golden_gate_assembly:v0"
    assert event.input_version_ids == tuple(fragment.version_id for fragment in fragments)
    assert event.output_version_id == product.version_id


def test_pcr_product_has_provenance_event() -> None:
    record = SequenceRecord(id="pcr", sequence="AAAACCCCGGGGTTTT")

    product = simulate_pcr(record, "AAAA", "AAAA", min_anneal_length=4)

    assert product.provenance is not None
    assert product.provenance.operation == "pcr"
    assert product.provenance.input_record_ids == ("pcr",)
    assert product.provenance.input_version_ids == (record.version_id,)
    assert product.provenance.output_record_id == "pcr_pcr_product"
    assert product.provenance.output_version_id is not None
    assert product.provenance.output_version_id.startswith("pcr_pcr_product:v")
    assert product.provenance.output_version_id != "pcr_pcr_product:v0"
    assert product.provenance.affected_ranges == ((0, 16),)
    assert product.provenance.parameters["length"] == product.length


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
