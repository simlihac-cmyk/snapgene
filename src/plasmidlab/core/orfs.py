"""ORF detection, translation, CDS validation, and reverse translation."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from Bio.Data import CodonTable

from plasmidlab.core.models import Feature, FeatureSegment, MoleculeTopology, MoleculeType, SequenceRecord


DNA_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


@dataclass(frozen=True, slots=True)
class ORF:
    """An open reading frame in original template coordinates."""

    start: int
    end: int
    strand: int
    frame: int
    segments: tuple[FeatureSegment, ...]
    nucleotide_sequence: str
    protein_sequence: str
    start_codon: str
    stop_codon: str
    wraps_origin: bool = False

    @property
    def nucleotide_length(self) -> int:
        """Return the ORF nucleotide length, including the stop codon."""

        return len(self.nucleotide_sequence)

    @property
    def amino_acid_length(self) -> int:
        """Return protein length excluding the terminal stop symbol."""

        return len(self.protein_sequence.rstrip("*"))


@dataclass(frozen=True, slots=True)
class CDSValidationResult:
    """Validation result for a candidate CDS."""

    is_valid: bool
    issues: tuple[str, ...]
    start_codon: str | None
    stop_codon: str | None
    internal_stop_positions: tuple[int, ...]
    nucleotide_length: int


def find_orfs(
    record: SequenceRecord,
    *,
    start_codons: tuple[str, ...] | None = ("ATG",),
    stop_codons: tuple[str, ...] | None = None,
    min_aa_length: int = 0,
    genetic_code: int | str = 1,
    include_reverse: bool = True,
) -> tuple[ORF, ...]:
    """Find closed ORFs on DNA in three or six frames."""

    _require_dna(record)
    if min_aa_length < 0:
        msg = "minimum amino acid length must be non-negative"
        raise ValueError(msg)
    table = _codon_table(genetic_code)
    starts = _normalize_codons(tuple(table.start_codons) if start_codons is None else start_codons)
    stops = _normalize_codons(tuple(table.stop_codons) if stop_codons is None else stop_codons)

    strands = (1, -1) if include_reverse else (1,)
    orfs: list[ORF] = []
    for strand in strands:
        oriented_sequence = record.sequence if strand == 1 else _reverse_complement(record.sequence)
        for frame in range(3):
            orfs.extend(
                _find_orfs_on_oriented_sequence(
                    record,
                    oriented_sequence,
                    strand=strand,
                    frame=frame,
                    start_codons=starts,
                    stop_codons=stops,
                    min_aa_length=min_aa_length,
                    genetic_code=genetic_code,
                )
            )
    return tuple(sorted(orfs, key=lambda orf: (orf.start, orf.strand, orf.nucleotide_length)))


def translate_sequence(
    sequence: str,
    *,
    genetic_code: int | str = 1,
    stop_symbol: str = "*",
) -> str:
    """Translate a DNA sequence string, ignoring an incomplete trailing codon."""

    normalized = _normalize_dna(sequence)
    table = _codon_table(genetic_code)
    amino_acids: list[str] = []
    for index in range(0, len(normalized) - len(normalized) % 3, 3):
        codon = normalized[index : index + 3]
        amino_acids.append(_translate_codon(codon, table, stop_symbol))
    return "".join(amino_acids)


def translate_record(
    record: SequenceRecord,
    *,
    start: int | None = None,
    end: int | None = None,
    strand: int = 1,
    genetic_code: int | str = 1,
    stop_symbol: str = "*",
) -> str:
    """Translate a whole record or selected region."""

    nucleotide_sequence = extract_region(record, start=start, end=end, strand=strand)
    return translate_sequence(nucleotide_sequence, genetic_code=genetic_code, stop_symbol=stop_symbol)


def translate_feature(
    record: SequenceRecord,
    feature: Feature,
    *,
    genetic_code: int | str = 1,
    stop_symbol: str = "*",
) -> str:
    """Translate a feature, including compound features."""

    nucleotide_sequence = _feature_sequence(record, feature)
    return translate_sequence(nucleotide_sequence, genetic_code=genetic_code, stop_symbol=stop_symbol)


def extract_region(
    record: SequenceRecord,
    *,
    start: int | None = None,
    end: int | None = None,
    strand: int = 1,
) -> str:
    """Extract a selected DNA region, with circular origin-crossing support."""

    _require_dna(record)
    if strand not in (-1, 1):
        msg = "translation strand must be +1 or -1"
        raise ValueError(msg)
    start = 0 if start is None else start
    end = record.length if end is None else end
    if start < 0 or end < 0 or start > record.length or end > record.length:
        msg = "region coordinates must be within the sequence"
        raise ValueError(msg)
    if record.topology is MoleculeTopology.LINEAR and start > end:
        msg = "linear regions cannot cross the origin"
        raise ValueError(msg)
    if start <= end:
        sequence = record.sequence[start:end]
    else:
        sequence = record.sequence[start:] + record.sequence[:end]
    return _reverse_complement(sequence) if strand == -1 else sequence


def validate_cds(
    record_or_sequence: SequenceRecord | str,
    *,
    feature: Feature | None = None,
    start: int | None = None,
    end: int | None = None,
    strand: int = 1,
    genetic_code: int | str = 1,
    start_codons: tuple[str, ...] | None = ("ATG",),
    stop_codons: tuple[str, ...] | None = None,
) -> CDSValidationResult:
    """Validate a candidate CDS for frame, start, stop, and internal stops."""

    if isinstance(record_or_sequence, SequenceRecord):
        if feature is not None:
            sequence = _feature_sequence(record_or_sequence, feature)
        else:
            sequence = extract_region(record_or_sequence, start=start, end=end, strand=strand)
    else:
        sequence = _normalize_dna(record_or_sequence)
        if strand == -1:
            sequence = _reverse_complement(sequence)
        elif strand != 1:
            msg = "CDS strand must be +1 or -1"
            raise ValueError(msg)

    table = _codon_table(genetic_code)
    starts = _normalize_codons(tuple(table.start_codons) if start_codons is None else start_codons)
    stops = _normalize_codons(tuple(table.stop_codons) if stop_codons is None else stop_codons)
    issues: list[str] = []
    start_codon = sequence[:3] if len(sequence) >= 3 else None
    stop_codon = sequence[-3:] if len(sequence) >= 3 else None

    if len(sequence) == 0:
        issues.append("CDS is empty")
    if len(sequence) % 3:
        issues.append("CDS length is not divisible by 3")
    if start_codon not in starts:
        issues.append("CDS does not start with an allowed start codon")
    if stop_codon not in stops:
        issues.append("CDS does not end with an allowed stop codon")

    internal_stops = tuple(
        index
        for index in range(3, max(len(sequence) - 3, 3), 3)
        if sequence[index : index + 3] in stops
    )
    if internal_stops:
        issues.append("CDS contains internal stop codons")

    return CDSValidationResult(
        is_valid=not issues,
        issues=tuple(issues),
        start_codon=start_codon,
        stop_codon=stop_codon,
        internal_stop_positions=internal_stops,
        nucleotide_length=len(sequence),
    )


def reverse_translate(
    protein_sequence: str,
    *,
    genetic_code: int | str = 1,
    codon_usage: Mapping[str, float] | None = None,
    codon_usage_path: str | Path | None = None,
) -> str:
    """Back-translate a protein sequence using a deterministic codon choice."""

    protein = _normalize_protein(protein_sequence)
    usage = dict(codon_usage or {})
    if codon_usage_path is not None:
        usage.update(load_codon_usage(codon_usage_path))

    table = _codon_table(genetic_code)
    codons_by_aa = _codons_by_amino_acid(table)
    dna: list[str] = []
    for amino_acid in protein:
        codons = codons_by_aa.get(amino_acid)
        if not codons:
            msg = f"cannot reverse translate unsupported amino acid: {amino_acid}"
            raise ValueError(msg)
        dna.append(_choose_codon(codons, usage))
    return "".join(dna)


def load_codon_usage(path: str | Path) -> dict[str, float]:
    """Load codon usage weights from JSON or CSV."""

    usage_path = Path(path)
    if usage_path.suffix.lower() == ".json":
        data = json.loads(usage_path.read_text(encoding="utf-8"))
        if isinstance(data, Mapping) and "codons" in data:
            data = data["codons"]
        if not isinstance(data, Mapping):
            msg = "JSON codon usage must be an object mapping codons to weights"
            raise ValueError(msg)
        return _normalize_usage_mapping(data)
    if usage_path.suffix.lower() == ".csv":
        with usage_path.open(newline="", encoding="utf-8") as handle:
            rows = csv.DictReader(handle)
            usage: dict[str, float] = {}
            for row in rows:
                codon = row.get("codon") or row.get("Codon")
                weight = (
                    row.get("weight")
                    or row.get("Weight")
                    or row.get("frequency")
                    or row.get("Frequency")
                    or row.get("usage")
                    or row.get("Usage")
                )
                if codon is None or weight is None:
                    msg = "CSV codon usage requires codon and weight/frequency columns"
                    raise ValueError(msg)
                usage[codon.upper().replace("U", "T")] = float(weight)
            return _normalize_usage_mapping(usage)
    msg = "codon usage path must end in .json or .csv"
    raise ValueError(msg)


def _find_orfs_on_oriented_sequence(
    record: SequenceRecord,
    oriented_sequence: str,
    *,
    strand: int,
    frame: int,
    start_codons: frozenset[str],
    stop_codons: frozenset[str],
    min_aa_length: int,
    genetic_code: int | str,
) -> list[ORF]:
    sequence_length = record.length
    if sequence_length < 6:
        return []
    search_sequence = (
        oriented_sequence + oriented_sequence
        if record.topology is MoleculeTopology.CIRCULAR
        else oriented_sequence
    )
    orfs: list[ORF] = []
    for position in range(frame, sequence_length, 3):
        codon = search_sequence[position : position + 3]
        if len(codon) < 3 or codon not in start_codons:
            continue
        scan_limit = (
            position + sequence_length
            if record.topology is MoleculeTopology.CIRCULAR
            else sequence_length
        )
        for stop_position in range(position + 3, scan_limit - 2, 3):
            stop_codon = search_sequence[stop_position : stop_position + 3]
            if stop_codon not in stop_codons:
                continue
            nucleotide_sequence = search_sequence[position : stop_position + 3]
            protein_sequence = translate_sequence(nucleotide_sequence, genetic_code=genetic_code)
            amino_acid_length = len(protein_sequence.rstrip("*"))
            if amino_acid_length >= min_aa_length:
                orfs.append(
                    _orf_from_oriented_interval(
                        record,
                        strand=strand,
                        frame=frame,
                        oriented_start=position,
                        oriented_end=stop_position + 3,
                        nucleotide_sequence=nucleotide_sequence,
                        protein_sequence=protein_sequence,
                        start_codon=codon,
                        stop_codon=stop_codon,
                    )
                )
            break
    return orfs


def _orf_from_oriented_interval(
    record: SequenceRecord,
    *,
    strand: int,
    frame: int,
    oriented_start: int,
    oriented_end: int,
    nucleotide_sequence: str,
    protein_sequence: str,
    start_codon: str,
    stop_codon: str,
) -> ORF:
    segments = _map_oriented_segments(record.length, strand, oriented_start, oriented_end)
    wraps_origin = any(segment.start == 0 for segment in segments[1:]) or (
        len(segments) > 1 and segments[0].start > segments[-1].start
    )
    return ORF(
        start=segments[0].start,
        end=segments[-1].end,
        strand=strand,
        frame=frame + 1 if strand == 1 else -(frame + 1),
        segments=segments,
        nucleotide_sequence=nucleotide_sequence,
        protein_sequence=protein_sequence,
        start_codon=start_codon,
        stop_codon=stop_codon,
        wraps_origin=wraps_origin,
    )


def _map_oriented_segments(
    sequence_length: int,
    strand: int,
    oriented_start: int,
    oriented_end: int,
) -> tuple[FeatureSegment, ...]:
    oriented_parts = _split_circular_interval(oriented_start, oriented_end, sequence_length)
    mapped: list[FeatureSegment] = []
    for part_start, part_end in oriented_parts:
        if strand == 1:
            mapped.append(FeatureSegment(part_start, part_end, strand=1))
        else:
            mapped.append(FeatureSegment(sequence_length - part_end, sequence_length - part_start, strand=-1))
    return tuple(mapped)


def _split_circular_interval(
    start: int,
    end: int,
    sequence_length: int,
) -> tuple[tuple[int, int], ...]:
    normalized_start = start % sequence_length
    if end <= sequence_length:
        return ((normalized_start, end),)
    normalized_end = end % sequence_length
    return ((normalized_start, sequence_length), (0, normalized_end))


def _feature_sequence(record: SequenceRecord, feature: Feature) -> str:
    _require_dna(record)
    parts = [record.sequence[segment.start : segment.end] for segment in feature.segments]
    sequence = "".join(parts)
    return _reverse_complement(sequence) if feature.strand == -1 else sequence


def _codon_table(genetic_code: int | str) -> Any:
    if isinstance(genetic_code, int):
        return CodonTable.unambiguous_dna_by_id[genetic_code]
    try:
        return CodonTable.unambiguous_dna_by_name[genetic_code]
    except KeyError:
        return CodonTable.unambiguous_dna_by_id[int(genetic_code)]


def _codons_by_amino_acid(table: Any) -> dict[str, tuple[str, ...]]:
    grouped: dict[str, list[str]] = defaultdict(list)
    for codon, amino_acid in table.forward_table.items():
        grouped[amino_acid].append(codon)
    grouped["*"].extend(table.stop_codons)
    return {amino_acid: tuple(sorted(codons)) for amino_acid, codons in grouped.items()}


def _choose_codon(codons: tuple[str, ...], usage: Mapping[str, float]) -> str:
    return sorted(codons, key=lambda codon: (-float(usage.get(codon, 0.0)), codon))[0]


def _translate_codon(codon: str, table: Any, stop_symbol: str) -> str:
    if codon in table.stop_codons:
        return stop_symbol
    amino_acid = table.forward_table.get(codon)
    if amino_acid is None:
        return "X"
    return amino_acid


def _normalize_codons(codons: tuple[str, ...]) -> frozenset[str]:
    normalized = frozenset(_normalize_dna(codon) for codon in codons)
    if any(len(codon) != 3 for codon in normalized):
        msg = "codons must be exactly three bases"
        raise ValueError(msg)
    return normalized


def _normalize_dna(sequence: str) -> str:
    if not isinstance(sequence, str):
        msg = "DNA sequence must be a string"
        raise TypeError(msg)
    normalized = sequence.upper().replace("U", "T")
    invalid = sorted(set(normalized) - set("ACGTRYSWKMBDHVN"))
    if invalid:
        msg = f"invalid DNA sequence characters: {''.join(invalid)}"
        raise ValueError(msg)
    return normalized


def _normalize_protein(sequence: str) -> str:
    if not isinstance(sequence, str):
        msg = "protein sequence must be a string"
        raise TypeError(msg)
    normalized = sequence.upper()
    invalid = sorted(set(normalized) - set("ACDEFGHIKLMNPQRSTVWY*"))
    if invalid:
        msg = f"invalid protein sequence characters: {''.join(invalid)}"
        raise ValueError(msg)
    return normalized


def _normalize_usage_mapping(values: Mapping[str, Any]) -> dict[str, float]:
    usage: dict[str, float] = {}
    for raw_codon, raw_weight in values.items():
        codon = _normalize_dna(str(raw_codon))
        if len(codon) != 3:
            msg = "codon usage keys must be codons"
            raise ValueError(msg)
        usage[codon] = float(raw_weight)
    return usage


def _reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(DNA_COMPLEMENT)[::-1]


def _require_dna(record: SequenceRecord) -> None:
    if record.molecule_type is not MoleculeType.DNA:
        msg = "ORF and translation tools require a DNA sequence"
        raise ValueError(msg)
