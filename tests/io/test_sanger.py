from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from plasmidlab.io.sanger import (
    SangerTrace,
    align_trace_to_reference,
    trace_from_abif_raw,
    trace_from_biopython_record,
)


def test_trace_from_mocked_abif_raw_extracts_calls_quality_peaks_and_channels() -> None:
    raw = {
        "PBAS2": b"ACGT",
        "PCON2": bytes([30, 31, 32, 33]),
        "PLOC2": [2, 5, 8, 11],
        "FWO_1": b"ACGT",
        "DATA9": [0, 20, 100, 20, 0],
        "DATA10": [0, 0, 10, 90, 10],
        "DATA11": [0, 0, 0, 10, 80],
        "DATA12": [70, 10, 0, 0, 0],
    }

    trace = trace_from_abif_raw(raw, id="mock")

    assert trace.id == "mock"
    assert trace.called_bases == "ACGT"
    assert trace.qualities == (30, 31, 32, 33)
    assert trace.peak_positions == (2, 5, 8, 11)
    assert trace.chromatogram["A"] == (0, 20, 100, 20, 0)
    assert trace.chromatogram["C"] == (0, 0, 10, 90, 10)
    assert trace.chromatogram["G"] == (0, 0, 0, 10, 80)
    assert trace.chromatogram["T"] == (70, 10, 0, 0, 0)


def test_trace_from_biopython_record_uses_sequence_and_letter_quality() -> None:
    bio_record = SeqRecord(Seq("ACGT"), id="abi_like")
    bio_record.letter_annotations["phred_quality"] = [35, 34, 33, 32]
    bio_record.annotations["abif_raw"] = {
        "PLOC2": [1, 4, 7, 10],
        "FWO_1": "GATC",
        "DATA9": [1, 2],
        "DATA10": [3, 4],
        "DATA11": [5, 6],
        "DATA12": [7, 8],
    }

    trace = trace_from_biopython_record(bio_record)

    assert trace.called_bases == "ACGT"
    assert trace.qualities == (35, 34, 33, 32)
    assert trace.peak_positions == (1, 4, 7, 10)
    assert trace.chromatogram["G"] == (1, 2)
    assert trace.chromatogram["A"] == (3, 4)
    assert trace.chromatogram["T"] == (5, 6)
    assert trace.chromatogram["C"] == (7, 8)


def test_align_trace_to_reference_reports_discrepancy_quality() -> None:
    trace = SangerTrace(id="trace", called_bases="ACCT", qualities=(40, 39, 20, 38))

    result = align_trace_to_reference("ACGT", trace)

    assert len(result.discrepancies) == 1
    assert result.discrepancies[0].position == 2
    assert result.discrepancies[0].reference_base == "G"
    assert result.discrepancies[0].observed_base == "C"
    assert result.discrepancies[0].quality == 20
