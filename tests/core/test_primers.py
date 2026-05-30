import pytest

from plasmidlab.core import (
    MoleculeTopology,
    PCRAmbiguousProductError,
    PCRSimulationError,
    SequenceRecord,
    design_primers,
    find_primer_bindings,
    primer_metrics,
    simulate_pcr,
)


def test_primer_binding_search_finds_exact_and_reverse_complement_sites() -> None:
    record = _circular_pcr_template()
    forward = record.sequence[60:72]
    reverse = _reverse_complement(record.sequence[10:22])

    forward_bindings = find_primer_bindings(record, forward)
    reverse_bindings = find_primer_bindings(record, reverse)

    assert [(binding.start, binding.end, binding.strand) for binding in forward_bindings] == [
        (60, 72, 1)
    ]
    assert [(binding.start, binding.end, binding.strand) for binding in reverse_bindings] == [
        (10, 22, -1)
    ]


def test_primer_binding_search_supports_mismatch_tolerance() -> None:
    record = SequenceRecord(id="mismatch", sequence="AACCGGTTAACC")

    assert find_primer_bindings(record, "AACCGGTTTACC") == ()

    bindings = find_primer_bindings(record, "AACCGGTTTACC", mismatches=1)

    assert len(bindings) == 1
    assert bindings[0].start == 0
    assert bindings[0].mismatches == 1


def test_primer_binding_search_supports_circular_origin_match() -> None:
    record = SequenceRecord(
        id="origin_binding",
        sequence="CCCGGG" + "A" * 20 + "AAATTT",
        topology=MoleculeTopology.CIRCULAR,
    )

    bindings = find_primer_bindings(record, "AAATTTCCCGGG")

    assert len(bindings) == 1
    assert bindings[0].start == 26
    assert bindings[0].end == 6
    assert bindings[0].wraps_origin


def test_primer_metrics_include_basic_and_primer3_values_when_available() -> None:
    metrics = primer_metrics("ATGCGCAT")

    assert metrics.length == 8
    assert metrics.gc_percent == 50.0
    assert metrics.wallace_tm == 24.0
    assert metrics.primer3_tm is None or isinstance(metrics.primer3_tm, float)
    assert isinstance(metrics.hairpin.available, bool)
    assert isinstance(metrics.self_dimer.available, bool)


def test_design_primers_wrapper_returns_primer3_pairs() -> None:
    pytest.importorskip("primer3")
    record = SequenceRecord(
        id="design",
        sequence=(
            "GCGTACGTAGCTAGCTAGCGTACGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
            "CGATCGATCGATCGATCGATCGATCGATCGATATGCGTACGTAGCTAGCTAGCGTACGATCGA"
            "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
        ),
    )

    pairs = design_primers(
        record,
        target_region=(40, 60),
        product_size_range=(80, 140),
        num_return=2,
    )

    assert pairs
    assert pairs[0].forward_sequence
    assert pairs[0].reverse_sequence
    assert 80 <= pairs[0].product_size <= 140


def test_pcr_across_origin_on_circular_plasmid() -> None:
    record = _circular_pcr_template()
    forward = record.sequence[60:72]
    reverse = _reverse_complement(record.sequence[10:22])

    product = simulate_pcr(record, forward, reverse, min_anneal_length=12)

    assert product.start == 60
    assert product.end == 22
    assert product.wraps_origin
    assert product.sequence == record.sequence[60:] + record.sequence[:22]
    assert product.length == 42


def test_pcr_includes_five_prime_primer_tails() -> None:
    record = _circular_pcr_template()
    forward_tail = "GGGG"
    reverse_tail = "CCCC"
    forward = forward_tail + record.sequence[60:72]
    reverse = reverse_tail + _reverse_complement(record.sequence[10:22])

    product = simulate_pcr(record, forward, reverse, min_anneal_length=12)

    assert product.sequence.startswith(forward_tail + record.sequence[60:72])
    assert product.sequence.endswith(record.sequence[10:22] + _reverse_complement(reverse_tail))
    assert product.length == 50
    assert product.forward_binding.tail_sequence == forward_tail
    assert product.reverse_binding.tail_sequence == reverse_tail
    assert product.provenance is not None
    assert product.provenance.parameters["forward_binding"]["anneal_length"] == 12
    assert product.provenance.parameters["forward_binding"]["tail_length"] == len(forward_tail)
    assert product.provenance.parameters["reverse_binding"]["tail_length"] == len(reverse_tail)


def test_tailed_primer_shorter_intended_match_is_not_blocked_by_longer_off_target() -> None:
    forward_tail = "GGGG"
    forward_anneal = "ATGACCTGATCG"
    reverse_site = "GCTAGTCCATGA"
    forward = forward_tail + forward_anneal
    reverse = _reverse_complement(reverse_site)
    off_target = forward + "A" * 8 + reverse_site
    spacer = "C" * 20
    intended_template = forward_anneal + "T" * 8 + reverse_site
    sequence = off_target + spacer + intended_template + "G" * 5
    intended_start = len(off_target) + len(spacer)
    intended_end = intended_start + len(intended_template)
    record = SequenceRecord(id="tail_off_target", sequence=sequence)

    with pytest.raises(PCRAmbiguousProductError, match="multiple possible PCR products"):
        simulate_pcr(record, forward, reverse, min_anneal_length=len(forward_anneal))

    product = simulate_pcr(
        record,
        forward,
        reverse,
        min_anneal_length=len(forward_anneal),
        max_anneal_length=len(forward),
        target_region=(intended_start, intended_end),
    )

    assert product.start == intended_start
    assert product.end == intended_end
    assert product.forward_binding.binding_length == len(forward_anneal)
    assert product.forward_binding.tail_sequence == forward_tail
    assert product.sequence == forward_tail + sequence[intended_start:intended_end]


def test_full_circle_circular_pcr_reports_explicit_nonzero_interval_and_tails() -> None:
    record = SequenceRecord(
        id="full_circle_pcr",
        sequence="AAAACCCCGGGGTTTT",
        topology=MoleculeTopology.CIRCULAR,
    )
    forward_tail = "GG"
    reverse_tail = "CC"
    forward = forward_tail + record.sequence[4:8]
    reverse = reverse_tail + _reverse_complement(record.sequence[0:4])

    product = simulate_pcr(record, forward, reverse, min_anneal_length=4)

    expected_template = record.sequence[4:] + record.sequence[:4]
    assert product.start == 4
    assert product.end == 4
    assert product.full_circle
    assert product.wraps_origin
    assert product.sequence == forward_tail + expected_template + _reverse_complement(reverse_tail)
    assert product.provenance is not None
    assert product.provenance.affected_ranges == ((4, record.length), (0, 4))
    assert product.provenance.parameters["full_circle"] is True
    assert product.provenance.parameters["template_interval_semantics"] == "full_circle"
    assert product.provenance.parameters["wraps_origin"] is True


def test_pcr_wrong_orientation_fails_clearly() -> None:
    record = _circular_pcr_template()
    forward = record.sequence[60:72]
    reverse = _reverse_complement(record.sequence[10:22])

    with pytest.raises(PCRSimulationError, match="do not face each other"):
        simulate_pcr(record, reverse, forward, min_anneal_length=12)


def test_pcr_multiple_binding_sites_are_ambiguous_unless_pair_selected() -> None:
    forward = "ATGACCTGATCG"
    reverse_site = "GCTAGTCCATGA"
    reverse = _reverse_complement(reverse_site)
    sequence = forward + "A" * 18 + forward + "C" * 18 + reverse_site + "G" * 10
    record = SequenceRecord(id="ambiguous", sequence=sequence)

    with pytest.raises(PCRAmbiguousProductError, match="multiple possible PCR products"):
        simulate_pcr(record, forward, reverse, min_anneal_length=12)

    selected = simulate_pcr(record, forward, reverse, min_anneal_length=12, pair_index=1)

    assert selected.start == 30
    assert selected.end == 72
    assert selected.sequence == sequence[30:72]


def _circular_pcr_template() -> SequenceRecord:
    sequence = (
        "TTGACCGTAA"
        "GCTAGTCCATGA"
        "CCGTAAGGTTCCGATACCGGTTACGATCGGATATCGAA"
        "ATGACCTGATCG"
        "GGTTAACC"
    )
    assert len(sequence) == 80
    return SequenceRecord(id="pcr_circle", sequence=sequence, topology=MoleculeTopology.CIRCULAR)


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
