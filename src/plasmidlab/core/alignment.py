"""Pairwise alignment and reference discrepancy reporting."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Sequence

from Bio import Align

from plasmidlab.core.models import MoleculeType, SequenceRecord


class AlignmentMode(StrEnum):
    """Supported pairwise alignment modes."""

    GLOBAL = "global"
    LOCAL = "local"
    SEMI_GLOBAL = "semi_global"


@dataclass(frozen=True, slots=True)
class PairwiseAlignmentResult:
    """A deterministic, serializable view of a pairwise alignment."""

    reference: str
    query: str
    aligned_reference: str
    aligned_query: str
    score: float
    mode: AlignmentMode
    molecule_type: MoleculeType
    reference_start: int
    reference_end: int
    query_start: int
    query_end: int


@dataclass(frozen=True, slots=True)
class Discrepancy:
    """A mismatch, insertion, or deletion relative to the reference sequence."""

    position: int
    kind: str
    reference_base: str
    observed_base: str
    quality: int | None = None


@dataclass(frozen=True, slots=True)
class ReferenceAlignmentResult:
    """Reference alignment plus discrepancy report."""

    alignment: PairwiseAlignmentResult
    discrepancies: tuple[Discrepancy, ...]


def pairwise_align(
    reference: SequenceRecord | str,
    query: SequenceRecord | str,
    *,
    mode: AlignmentMode | str = AlignmentMode.GLOBAL,
    molecule_type: MoleculeType | str | None = None,
    match_score: float = 2.0,
    mismatch_score: float = -1.0,
    open_gap_score: float = -2.0,
    extend_gap_score: float = -0.5,
) -> PairwiseAlignmentResult:
    """Align ``query`` to ``reference`` and return gapped alignment strings."""

    resolved_mode = _coerce_mode(mode)
    resolved_molecule_type = _coerce_molecule_type(molecule_type, reference, query)
    reference_sequence = _sequence_text(reference)
    query_sequence = _sequence_text(query)

    aligner = Align.PairwiseAligner()
    aligner.mode = "local" if resolved_mode is AlignmentMode.LOCAL else "global"
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score
    if resolved_mode is AlignmentMode.SEMI_GLOBAL:
        # Free terminal gaps make this a practical end-gap-free reference workflow.
        aligner.end_insertion_score = 0.0
        aligner.end_deletion_score = 0.0

    alignments = aligner.align(reference_sequence, query_sequence)
    if not alignments:
        msg = "no alignment was produced"
        raise ValueError(msg)
    alignment = alignments[0]
    aligned_reference, aligned_query = _gapped_strings(reference_sequence, query_sequence, alignment)
    coordinates = alignment.coordinates
    return PairwiseAlignmentResult(
        reference=reference_sequence,
        query=query_sequence,
        aligned_reference=aligned_reference,
        aligned_query=aligned_query,
        score=float(alignment.score),
        mode=resolved_mode,
        molecule_type=resolved_molecule_type,
        reference_start=int(coordinates[0][0]),
        reference_end=int(coordinates[0][-1]),
        query_start=int(coordinates[1][0]),
        query_end=int(coordinates[1][-1]),
    )


def align_to_reference(
    reference: SequenceRecord | str,
    observed: SequenceRecord | str,
    *,
    qualities: Sequence[int] | None = None,
    mode: AlignmentMode | str = AlignmentMode.GLOBAL,
    molecule_type: MoleculeType | str | None = None,
) -> ReferenceAlignmentResult:
    """Align an observed sequence to a reference and report base-level differences."""

    alignment = pairwise_align(
        reference,
        observed,
        mode=mode,
        molecule_type=molecule_type,
    )
    return ReferenceAlignmentResult(
        alignment=alignment,
        discrepancies=_discrepancies_from_alignment(alignment, qualities),
    )


def discrepancy_report_tsv(result: ReferenceAlignmentResult) -> str:
    """Return a tab-separated discrepancy report."""

    lines = ["position\ttype\treference\tobserved\tquality"]
    for discrepancy in result.discrepancies:
        quality = "" if discrepancy.quality is None else str(discrepancy.quality)
        lines.append(
            "\t".join(
                (
                    str(discrepancy.position),
                    discrepancy.kind,
                    discrepancy.reference_base,
                    discrepancy.observed_base,
                    quality,
                )
            )
        )
    return "\n".join(lines)


def _discrepancies_from_alignment(
    alignment: PairwiseAlignmentResult,
    qualities: Sequence[int] | None,
) -> tuple[Discrepancy, ...]:
    discrepancies: list[Discrepancy] = []
    reference_position = alignment.reference_start
    query_position = alignment.query_start
    for reference_base, observed_base in zip(
        alignment.aligned_reference,
        alignment.aligned_query,
        strict=True,
    ):
        quality = _quality_at(qualities, query_position) if observed_base != "-" else None
        if reference_base == "-" and observed_base != "-":
            discrepancies.append(
                Discrepancy(
                    position=reference_position,
                    kind="insertion",
                    reference_base="-",
                    observed_base=observed_base,
                    quality=quality,
                )
            )
            query_position += 1
            continue
        if reference_base != "-" and observed_base == "-":
            discrepancies.append(
                Discrepancy(
                    position=reference_position,
                    kind="deletion",
                    reference_base=reference_base,
                    observed_base="-",
                    quality=None,
                )
            )
            reference_position += 1
            continue
        if reference_base.upper() != observed_base.upper():
            discrepancies.append(
                Discrepancy(
                    position=reference_position,
                    kind="mismatch",
                    reference_base=reference_base,
                    observed_base=observed_base,
                    quality=quality,
                )
            )
        reference_position += 1
        query_position += 1
    return tuple(discrepancies)


def _gapped_strings(reference: str, query: str, alignment: object) -> tuple[str, str]:
    coordinates = alignment.coordinates
    aligned_reference: list[str] = []
    aligned_query: list[str] = []
    for index in range(coordinates.shape[1] - 1):
        reference_start = int(coordinates[0][index])
        reference_end = int(coordinates[0][index + 1])
        query_start = int(coordinates[1][index])
        query_end = int(coordinates[1][index + 1])
        reference_length = reference_end - reference_start
        query_length = query_end - query_start
        if reference_length and query_length:
            aligned_reference.append(reference[reference_start:reference_end])
            aligned_query.append(query[query_start:query_end])
        elif reference_length:
            aligned_reference.append(reference[reference_start:reference_end])
            aligned_query.append("-" * reference_length)
        elif query_length:
            aligned_reference.append("-" * query_length)
            aligned_query.append(query[query_start:query_end])
    return "".join(aligned_reference), "".join(aligned_query)


def _quality_at(qualities: Sequence[int] | None, query_position: int) -> int | None:
    if qualities is None or query_position < 0 or query_position >= len(qualities):
        return None
    return int(qualities[query_position])


def _sequence_text(record_or_sequence: SequenceRecord | str) -> str:
    if isinstance(record_or_sequence, SequenceRecord):
        return record_or_sequence.sequence
    return "".join(str(record_or_sequence).split()).upper()


def _coerce_mode(mode: AlignmentMode | str) -> AlignmentMode:
    try:
        return mode if isinstance(mode, AlignmentMode) else AlignmentMode(mode)
    except ValueError as error:
        msg = f"unsupported alignment mode: {mode!r}"
        raise ValueError(msg) from error


def _coerce_molecule_type(
    molecule_type: MoleculeType | str | None,
    reference: SequenceRecord | str,
    query: SequenceRecord | str,
) -> MoleculeType:
    if molecule_type is not None:
        return molecule_type if isinstance(molecule_type, MoleculeType) else MoleculeType(molecule_type)
    if isinstance(reference, SequenceRecord):
        return reference.molecule_type
    if isinstance(query, SequenceRecord):
        return query.molecule_type
    return MoleculeType.DNA
