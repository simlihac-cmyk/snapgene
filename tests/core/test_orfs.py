import json

from plasmidlab.core import (
    Feature,
    MoleculeTopology,
    SequenceRecord,
    extract_region,
    find_orfs,
    load_codon_usage,
    reverse_translate,
    translate_feature,
    translate_record,
    translate_sequence,
    validate_cds,
)


def test_find_orfs_detects_plus_strand_orf() -> None:
    record = SequenceRecord(id="plus_orf", sequence="CCCATGAAATAAGGG")

    orfs = find_orfs(record, min_aa_length=2)

    assert [(orf.start, orf.end, orf.strand, orf.frame, orf.protein_sequence) for orf in orfs] == [
        (3, 12, 1, 1, "MK*")
    ]


def test_find_orfs_detects_minus_strand_orf() -> None:
    record = SequenceRecord(id="minus_orf", sequence="CCCTTATTTCATGGG")

    orfs = find_orfs(record, min_aa_length=2)

    minus_orfs = [orf for orf in orfs if orf.strand == -1]
    assert len(minus_orfs) == 1
    assert minus_orfs[0].start == 3
    assert minus_orfs[0].end == 12
    assert minus_orfs[0].frame == -1
    assert minus_orfs[0].nucleotide_sequence == "ATGAAATAA"
    assert minus_orfs[0].protein_sequence == "MK*"


def test_find_orfs_detects_circular_orf_crossing_origin() -> None:
    record = SequenceRecord(
        id="circular_orf",
        sequence="AATAA" + "C" * 9 + "ATGA",
        topology=MoleculeTopology.CIRCULAR,
    )

    orfs = find_orfs(record, min_aa_length=2, include_reverse=False)

    assert len(orfs) == 1
    assert orfs[0].start == 14
    assert orfs[0].end == 5
    assert orfs[0].wraps_origin
    assert [(segment.start, segment.end, segment.strand) for segment in orfs[0].segments] == [
        (14, 18, 1),
        (0, 5, 1),
    ]
    assert orfs[0].nucleotide_sequence == "ATGAAATAA"


def test_find_orfs_supports_configurable_start_and_stop_codons() -> None:
    record = SequenceRecord(id="alt_start", sequence="CCCGTGAAATAAGGG")

    default_orfs = find_orfs(record, min_aa_length=2, include_reverse=False)
    custom_orfs = find_orfs(
        record,
        start_codons=("GTG",),
        stop_codons=("TAA",),
        min_aa_length=2,
        include_reverse=False,
    )

    assert default_orfs == ()
    assert custom_orfs[0].start_codon == "GTG"
    assert custom_orfs[0].protein_sequence == "VK*"


def test_translation_supports_sequence_record_region_feature_and_code_table() -> None:
    record = SequenceRecord(id="translate", sequence="CCCATGAAATAAGGG")
    feature = Feature(type="CDS", start=3, end=12, strand=1, qualifiers={"gene": "mk"})
    reverse_record = SequenceRecord(id="reverse_translate", sequence="CCCTTATTTCATGGG")
    reverse_feature = Feature(type="CDS", start=3, end=12, strand=-1)

    assert translate_sequence("ATGAAATAA") == "MK*"
    assert translate_record(record, start=3, end=12) == "MK*"
    assert translate_record(SequenceRecord(id="whole", sequence="ATGTGA"), genetic_code=2) == "MW"
    assert translate_feature(record, feature) == "MK*"
    assert translate_feature(reverse_record, reverse_feature) == "MK*"
    assert extract_region(record, start=3, end=12) == "ATGAAATAA"


def test_cds_validation_reports_start_stop_internal_stop_and_frame() -> None:
    valid = validate_cds("ATGAAATAA")
    internal_stop = validate_cds("ATGTAATAA")
    bad_frame = validate_cds("ATGAAAAC")
    missing_stop = validate_cds("ATGAAA")

    assert valid.is_valid
    assert internal_stop.issues == ("CDS contains internal stop codons",)
    assert internal_stop.internal_stop_positions == (3,)
    assert "CDS length is not divisible by 3" in bad_frame.issues
    assert "CDS does not end with an allowed stop codon" in missing_stop.issues


def test_reverse_translation_uses_table_and_codon_usage() -> None:
    assert reverse_translate("MA") == "ATGGCA"
    assert reverse_translate("M*", codon_usage={"TGA": 10.0, "TAA": 1.0}) == "ATGTGA"


def test_load_codon_usage_from_json_and_csv(tmp_path) -> None:
    json_path = tmp_path / "usage.json"
    csv_path = tmp_path / "usage.csv"
    json_path.write_text(json.dumps({"codons": {"GCT": 1.0, "GCC": 5.0}}), encoding="utf-8")
    csv_path.write_text("codon,frequency\nTAA,1\nTAG,3\n", encoding="utf-8")

    assert load_codon_usage(json_path) == {"GCT": 1.0, "GCC": 5.0}
    assert load_codon_usage(csv_path) == {"TAA": 1.0, "TAG": 3.0}
    assert reverse_translate("A", codon_usage_path=json_path) == "GCC"
