from Bio.Restriction import BamHI, EcoRI

from plasmidlab.core import (
    FragmentEndSide,
    MoleculeTopology,
    Overhang,
    OverhangKind,
    RestrictionEnzyme,
    SequenceRecord,
    analyze_restriction_sites,
    compatible_fragment_ends,
    digest,
    fragment_end_for_site,
    find_restriction_sites,
)


def test_enzyme_search_classifies_non_single_and_multi_cutters() -> None:
    record = _engineered_plasmid(sequence_suffix="GAATTC")

    analysis = analyze_restriction_sites(record, "EcoRI,BamHI,HindIII,SmaI")

    assert analysis.non_cutters == ("SmaI",)
    assert set(analysis.single_cutters) == {"BamHI", "HindIII"}
    assert analysis.multi_cutters == ("EcoRI",)
    assert [site.start for site in analysis.sites_by_enzyme["BamHI"]] == [40]
    assert analysis.sites_by_enzyme["EcoRI"][0].recognition_sequence == "GAATTC"
    assert analysis.sites_by_enzyme["EcoRI"][0].overhang.kind is OverhangKind.FIVE_PRIME
    assert analysis.sites_by_enzyme["EcoRI"][0].overhang.sequence == "AATT"
    assert analysis.sites_by_enzyme["EcoRI"][0].methylation_sensitive is None


def test_circular_double_digest_returns_expected_fragment_sizes() -> None:
    record = _engineered_plasmid()

    fragments = digest(record, "EcoRI,BamHI")

    assert [(fragment.start, fragment.end, fragment.length) for fragment in fragments] == [
        (11, 41, 30),
        (41, 11, 70),
    ]
    assert [fragment.sequence for fragment in fragments] == [
        record.sequence[11:41],
        record.sequence[41:] + record.sequence[:11],
    ]
    assert fragments[0].source_enzymes == ("EcoRI", "BamHI")
    assert sum(fragment.length for fragment in fragments) == record.length


def test_circular_single_digest_returns_one_linearized_fragment() -> None:
    record = _engineered_plasmid()

    fragments = digest(record, "HindIII")

    assert len(fragments) == 1
    assert fragments[0].start == 71
    assert fragments[0].end == 71
    assert fragments[0].length == record.length
    assert fragments[0].sequence == record.sequence[71:] + record.sequence[:71]


def test_linear_digest_returns_edge_fragments() -> None:
    record = SequenceRecord(
        id="linear_digest",
        sequence=_engineered_sequence(),
        topology=MoleculeTopology.LINEAR,
    )

    fragments = digest(record, "EcoRI")

    assert [(fragment.start, fragment.end, fragment.length) for fragment in fragments] == [
        (0, 11, 11),
        (11, 100, 89),
    ]
    assert fragments[0].left_overhang.kind is OverhangKind.NONE
    assert fragments[0].right_overhang.kind is OverhangKind.FIVE_PRIME
    assert fragments[1].left_overhang.sequence == "AATT"
    assert fragments[1].right_overhang.kind is OverhangKind.NONE
    assert fragments[0].left_end.source_enzyme is None
    assert fragments[1].right_end.source_enzyme is None


def test_multi_enzyme_digest_accepts_biopython_classes_and_custom_sets() -> None:
    record = _engineered_plasmid()

    class_fragments = digest(record, (EcoRI, BamHI))
    named_fragments = digest(record, "EcoRI,BamHI")
    custom = RestrictionEnzyme(
        name="CustomEco",
        recognition_sequence="GAATTC",
        top_cut_offset=1,
        bottom_cut_offset=5,
        overhang=Overhang(OverhangKind.FIVE_PRIME, "AATT"),
        methylation_sensitive=True,
    )
    custom_sites = find_restriction_sites(record, (custom,))

    assert [fragment.length for fragment in class_fragments] == [
        fragment.length for fragment in named_fragments
    ]
    assert len(custom_sites) == 1
    assert custom_sites[0].enzyme_name == "CustomEco"
    assert custom_sites[0].methylation_sensitive is True


def test_enzyme_search_detects_reverse_complement_recognition_site() -> None:
    sequence = "A" * 20 + "GAGACC" + "C" * 20
    record = SequenceRecord(id="reverse_site", sequence=sequence)

    sites = find_restriction_sites(record, "BsaI")

    assert len(sites) == 1
    assert sites[0].start == 20
    assert sites[0].end == 26
    assert sites[0].strand == -1
    assert sites[0].top_cut == 15
    assert sites[0].bottom_cut == 19


def test_bsa_i_type_iis_digest_uses_sequence_derived_overhang_not_nnnn() -> None:
    record = SequenceRecord(
        id="bsa_i_overhang",
        sequence="TTTTGGTCTCAATGCGGGG",
        topology=MoleculeTopology.LINEAR,
    )

    fragments = digest(record, "BsaI")

    assert [(fragment.start, fragment.end) for fragment in fragments] == [(0, 11), (11, 19)]
    assert fragments[0].right_end.kind is OverhangKind.FIVE_PRIME
    assert fragments[0].right_end.top_strand_overhang_sequence == "ATGC"
    assert fragments[0].right_end.overhang_sequence == "GCAT"
    assert fragments[1].left_end.overhang_sequence == "ATGC"
    assert fragments[0].right_end.overhang_sequence != "NNNN"
    assert compatible_fragment_ends(fragments[0].right_end, fragments[1].left_end)


def test_bsa_i_derived_non_complementary_overhangs_do_not_ligate() -> None:
    first = digest(
        SequenceRecord(id="bsa_a", sequence="TTTTGGTCTCAATGCGGGG"),
        "BsaI",
    )
    second = digest(
        SequenceRecord(id="bsa_b", sequence="TTTTGGTCTCATTTTGGGG"),
        "BsaI",
    )

    assert not compatible_fragment_ends(first[0].right_end, second[1].left_end)


def test_ambiguous_overhangs_require_explicit_ambiguous_matching() -> None:
    ambiguous = digest(
        SequenceRecord(id="bsa_n", sequence="TTTTGGTCTCAANNNGGGG"),
        "BsaI",
    )
    concrete = digest(
        SequenceRecord(id="bsa_c", sequence="TTTTGGTCTCAATGCGGGG"),
        "BsaI",
    )

    assert not compatible_fragment_ends(ambiguous[0].right_end, concrete[1].left_end)
    assert compatible_fragment_ends(
        ambiguous[0].right_end,
        concrete[1].left_end,
        allow_ambiguous=True,
    )


def test_reversed_orientation_end_does_not_pass_as_compatible() -> None:
    fragments = digest(
        SequenceRecord(id="bsa_i_orientation", sequence="TTTTGGTCTCAATGCGGGG"),
        "BsaI",
    )

    assert not compatible_fragment_ends(fragments[1].left_end, fragments[1].left_end)
    assert not compatible_fragment_ends(fragments[0].right_end, fragments[0].right_end)


def test_classic_sticky_ends_are_derived_from_template_sequence() -> None:
    cases = (
        ("EcoRI", "AAAGAATTCTTT", "AATT"),
        ("BamHI", "AAAGGATCCTTT", "GATC"),
        ("HindIII", "AAAAAGCTTTTT", "AGCT"),
    )

    for enzyme, sequence, expected in cases:
        fragments = digest(SequenceRecord(id=enzyme, sequence=sequence), enzyme)
        assert fragments[0].right_end.kind is OverhangKind.FIVE_PRIME
        assert fragments[0].right_end.top_strand_overhang_sequence == expected
        assert fragments[1].left_end.top_strand_overhang_sequence == expected
        assert compatible_fragment_ends(fragments[0].right_end, fragments[1].left_end)


def test_blunt_end_compatibility_only_accepts_blunt_ends() -> None:
    blunt = digest(SequenceRecord(id="sma_i", sequence="AAACCCGGGTTT"), "SmaI")
    sticky = digest(SequenceRecord(id="eco_ri", sequence="AAAGAATTCTTT"), "EcoRI")

    assert blunt[0].right_end.kind is OverhangKind.BLUNT
    assert blunt[1].left_end.kind is OverhangKind.BLUNT
    assert compatible_fragment_ends(blunt[0].right_end, blunt[1].left_end)
    assert not compatible_fragment_ends(blunt[0].right_end, sticky[1].left_end)


def test_cut_geometry_records_actual_top_and_bottom_nicks() -> None:
    record = SequenceRecord(id="geometry", sequence="TTTTGGTCTCAATGCGGGG")
    site = find_restriction_sites(record, "BsaI")[0]

    geometry = site.cut_geometry
    left_end = fragment_end_for_site(record, site, FragmentEndSide.LEFT)

    assert geometry.top_cut == 11
    assert geometry.bottom_cut == 15
    assert geometry.cut_source_metadata["top_cut_absolute"] == 11
    assert left_end.cut_geometry == geometry


def test_circular_search_detects_origin_spanning_site() -> None:
    record = SequenceRecord(
        id="origin_site",
        sequence="TTC" + "A" * 20 + "GAA",
        topology=MoleculeTopology.CIRCULAR,
    )

    sites = find_restriction_sites(record, "EcoRI")

    assert len(sites) == 1
    assert sites[0].start == 23
    assert sites[0].end == 3
    assert sites[0].wraps_origin
    assert sites[0].top_cut == 24
    assert sites[0].bottom_cut == 2


def _engineered_plasmid(*, sequence_suffix: str = "") -> SequenceRecord:
    return SequenceRecord(
        id="pRE_digest",
        sequence=_engineered_sequence(sequence_suffix=sequence_suffix),
        topology=MoleculeTopology.CIRCULAR,
    )


def _engineered_sequence(*, sequence_suffix: str = "") -> str:
    base = (
        "A" * 10
        + "GAATTC"
        + "C" * 24
        + "GGATCC"
        + "T" * 24
        + "AAGCTT"
        + "G" * 24
    )
    if not sequence_suffix:
        return base
    return base + sequence_suffix
