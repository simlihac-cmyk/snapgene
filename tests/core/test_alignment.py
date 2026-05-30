from plasmidlab.core import (
    AlignmentMode,
    MoleculeType,
    SequenceRecord,
    align_to_reference,
    discrepancy_report_tsv,
    pairwise_align,
)
from plasmidlab.io.fasta import loads_fasta
from plasmidlab.io.genbank import dumps_genbank, loads_genbank


def test_global_alignment_reports_mismatch_with_quality() -> None:
    result = align_to_reference("ACGT", "AGGT", qualities=[30, 31, 32, 33])

    assert result.alignment.aligned_reference == "ACGT"
    assert result.alignment.aligned_query == "AGGT"
    assert result.discrepancies[0].kind == "mismatch"
    assert result.discrepancies[0].position == 1
    assert result.discrepancies[0].reference_base == "C"
    assert result.discrepancies[0].observed_base == "G"
    assert result.discrepancies[0].quality == 31


def test_global_alignment_reports_insertions_and_deletions() -> None:
    insertion = align_to_reference("ACGT", "ACGTA", qualities=[30, 31, 32, 33, 34])
    deletion = align_to_reference("ACGT", "AGT", qualities=[30, 31, 32])

    assert insertion.discrepancies[-1].kind == "insertion"
    assert insertion.discrepancies[-1].position == 4
    assert insertion.discrepancies[-1].reference_base == "-"
    assert insertion.discrepancies[-1].observed_base == "A"
    assert insertion.discrepancies[-1].quality == 34
    assert deletion.discrepancies[0].kind == "deletion"
    assert deletion.discrepancies[0].position == 1
    assert deletion.discrepancies[0].reference_base == "C"
    assert deletion.discrepancies[0].observed_base == "-"
    assert deletion.discrepancies[0].quality is None


def test_local_alignment_supports_protein_sequences() -> None:
    result = pairwise_align(
        "MEEPQSDPSV",
        "PQSD",
        mode=AlignmentMode.LOCAL,
        molecule_type=MoleculeType.PROTEIN,
    )

    assert result.mode is AlignmentMode.LOCAL
    assert result.molecule_type is MoleculeType.PROTEIN
    assert result.reference_start == 3
    assert result.reference_end == 7
    assert result.aligned_reference == "PQSD"
    assert result.aligned_query == "PQSD"


def test_semi_global_alignment_ignores_terminal_gaps() -> None:
    result = pairwise_align("AAAACCCCGGGG", "CCCC", mode=AlignmentMode.SEMI_GLOBAL)

    assert result.mode is AlignmentMode.SEMI_GLOBAL
    assert result.score == 8.0
    assert result.aligned_reference == "AAAACCCCGGGG"
    assert result.aligned_query == "----CCCC----"


def test_reference_alignment_accepts_fasta_and_genbank_records() -> None:
    fasta_reference = loads_fasta(">ref\nACGTACGT\n")[0]
    genbank_observed = loads_genbank(
        dumps_genbank(SequenceRecord(id="obs", sequence="ACGTTCGT"))
    )[0]

    result = align_to_reference(fasta_reference, genbank_observed)

    assert result.discrepancies[0].kind == "mismatch"
    assert result.discrepancies[0].position == 4
    assert result.discrepancies[0].reference_base == "A"
    assert result.discrepancies[0].observed_base == "T"


def test_discrepancy_report_tsv_lists_required_fields() -> None:
    result = align_to_reference("ACGT", "AGGT", qualities=[30, 31, 32, 33])

    report = discrepancy_report_tsv(result)

    assert report.splitlines()[0] == "position\ttype\treference\tobserved\tquality"
    assert "1\tmismatch\tC\tG\t31" in report
