"""Sanger AB1 trace import and reference alignment helpers."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord as BioSeqRecord

from plasmidlab.core.alignment import AlignmentMode, ReferenceAlignmentResult, align_to_reference
from plasmidlab.core.models import SequenceRecord


TRACE_BASES = ("A", "C", "G", "T")
ABI_TRACE_DATA_KEYS = ("DATA9", "DATA10", "DATA11", "DATA12")


@dataclass(frozen=True, slots=True)
class SangerTrace:
    """Called bases, qualities, and chromatogram signal from an AB1 trace."""

    id: str
    called_bases: str
    qualities: tuple[int, ...] = ()
    peak_positions: tuple[int, ...] = ()
    chromatogram: Mapping[str, tuple[int, ...]] = field(default_factory=dict)
    source_path: str | None = None

    def __post_init__(self) -> None:
        called_bases = "".join(self.called_bases.split()).upper()
        invalid = sorted(set(called_bases) - set("ACGTRYSWKMBDHVN"))
        if invalid:
            msg = f"invalid called bases in trace: {''.join(invalid)}"
            raise ValueError(msg)
        object.__setattr__(self, "called_bases", called_bases)
        object.__setattr__(self, "qualities", tuple(int(value) for value in self.qualities))
        object.__setattr__(self, "peak_positions", tuple(int(value) for value in self.peak_positions))
        object.__setattr__(
            self,
            "chromatogram",
            MappingProxyType(
                {
                    base.upper(): tuple(int(value) for value in values)
                    for base, values in self.chromatogram.items()
                }
            ),
        )


def read_ab1(path: str | Path) -> SangerTrace:
    """Read a Sanger AB1 file using Biopython's ABI parser."""

    source = Path(path)
    bio_record = SeqIO.read(source, "abi")
    trace = trace_from_biopython_record(bio_record)
    return SangerTrace(
        id=trace.id,
        called_bases=trace.called_bases,
        qualities=trace.qualities,
        peak_positions=trace.peak_positions,
        chromatogram=trace.chromatogram,
        source_path=str(source),
    )


def trace_from_biopython_record(bio_record: BioSeqRecord) -> SangerTrace:
    """Build a trace from a Biopython ABI SeqRecord."""

    abif_raw = bio_record.annotations.get("abif_raw", {})
    called_bases = str(bio_record.seq).upper()
    if not called_bases and isinstance(abif_raw, Mapping):
        called_bases = _decode_called_bases(abif_raw.get("PBAS2") or abif_raw.get("PBAS1"))
    qualities = tuple(int(value) for value in bio_record.letter_annotations.get("phred_quality", ()))
    if not qualities and isinstance(abif_raw, Mapping):
        qualities = _decode_quality_values(abif_raw.get("PCON2") or abif_raw.get("PCON1"))
    peak_positions: tuple[int, ...] = ()
    chromatogram: dict[str, tuple[int, ...]] = {}
    if isinstance(abif_raw, Mapping):
        peak_positions = _decode_ints(abif_raw.get("PLOC2") or abif_raw.get("PLOC1"))
        chromatogram = _decode_chromatogram(abif_raw)
    return SangerTrace(
        id=bio_record.id or bio_record.name or "trace",
        called_bases=called_bases,
        qualities=qualities,
        peak_positions=peak_positions,
        chromatogram=chromatogram,
    )


def trace_from_abif_raw(
    abif_raw: Mapping[str, object],
    *,
    id: str = "trace",
    source_path: str | None = None,
) -> SangerTrace:
    """Build a trace from an ABIF raw dictionary, useful for parser unit tests."""

    return SangerTrace(
        id=id,
        called_bases=_decode_called_bases(abif_raw.get("PBAS2") or abif_raw.get("PBAS1")),
        qualities=_decode_quality_values(abif_raw.get("PCON2") or abif_raw.get("PCON1")),
        peak_positions=_decode_ints(abif_raw.get("PLOC2") or abif_raw.get("PLOC1")),
        chromatogram=_decode_chromatogram(abif_raw),
        source_path=source_path,
    )


def align_trace_to_reference(
    reference: SequenceRecord | str,
    trace: SangerTrace,
    *,
    mode: AlignmentMode | str = AlignmentMode.GLOBAL,
) -> ReferenceAlignmentResult:
    """Align AB1 called bases to a reference and include trace quality scores."""

    return align_to_reference(
        reference,
        trace.called_bases,
        qualities=trace.qualities or None,
        mode=mode,
    )


def _decode_called_bases(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("ascii", errors="ignore").replace("\x00", "").strip().upper()
    if isinstance(value, str):
        return value.replace("\x00", "").strip().upper()
    if isinstance(value, Sequence):
        characters: list[str] = []
        for item in value:
            if isinstance(item, int):
                characters.append(chr(item))
            else:
                characters.append(str(item))
        return "".join(characters).replace("\x00", "").strip().upper()
    return str(value).replace("\x00", "").strip().upper()


def _decode_quality_values(value: object) -> tuple[int, ...]:
    if value is None:
        return ()
    if isinstance(value, bytes):
        return tuple(int(byte) for byte in value)
    return _decode_ints(value)


def _decode_ints(value: object) -> tuple[int, ...]:
    if value is None:
        return ()
    if isinstance(value, bytes):
        return tuple(int(byte) for byte in value)
    if isinstance(value, Sequence) and not isinstance(value, str):
        return tuple(int(item) for item in value)
    return (int(value),)


def _decode_chromatogram(abif_raw: Mapping[str, object]) -> dict[str, tuple[int, ...]]:
    channel_order = _channel_order(abif_raw.get("FWO_1"))
    chromatogram: dict[str, tuple[int, ...]] = {}
    for base, data_key in zip(channel_order, ABI_TRACE_DATA_KEYS, strict=True):
        values = _decode_ints(abif_raw.get(data_key))
        if values:
            chromatogram[base] = values
    return chromatogram


def _channel_order(value: object) -> tuple[str, str, str, str]:
    order = _decode_called_bases(value) if value is not None else "GATC"
    order = "".join(base for base in order if base in TRACE_BASES)
    if len(order) != 4 or set(order) != set(TRACE_BASES):
        return ("G", "A", "T", "C")
    return tuple(order)  # type: ignore[return-value]
