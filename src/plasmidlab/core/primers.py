"""Primer binding, metrics, design, and PCR simulation."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any

from plasmidlab.core.models import (
    MoleculeTopology,
    MoleculeType,
    ProvenanceEvent,
    SequenceRecord,
    _derive_version_id,
)

try:  # pragma: no cover - exercised when primer3-py is absent.
    import primer3.bindings as primer3_bindings
except ImportError:  # pragma: no cover
    primer3_bindings = None


IUPAC_DNA: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("GC"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}
DNA_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


class PrimerError(ValueError):
    """Base exception for primer analysis failures."""


class Primer3UnavailableError(PrimerError):
    """Raised when primer3-py is required but unavailable."""


class PCRSimulationError(PrimerError):
    """Raised when PCR simulation cannot produce a clear product."""


class PCRAmbiguousProductError(PCRSimulationError):
    """Raised when a primer pair can produce multiple products."""


@dataclass(frozen=True, slots=True)
class SecondaryStructureCheck:
    """Result from a primer3 secondary-structure check."""

    available: bool
    structure_found: bool | None = None
    tm: float | None = None
    dg: float | None = None


@dataclass(frozen=True, slots=True)
class PrimerMetrics:
    """Basic and primer3-backed primer metrics."""

    sequence: str
    length: int
    gc_percent: float
    wallace_tm: float
    primer3_tm: float | None
    hairpin: SecondaryStructureCheck
    self_dimer: SecondaryStructureCheck


@dataclass(frozen=True, slots=True)
class PrimerBinding:
    """A primer binding site on a template."""

    primer_sequence: str
    binding_sequence: str
    tail_sequence: str
    start: int
    end: int
    strand: int
    mismatches: int = 0
    wraps_origin: bool = False
    annealing_tm: float | None = None
    is_unique: bool = True

    def __post_init__(self) -> None:
        if self.strand not in (-1, 1):
            msg = "primer binding strand must be +1 or -1"
            raise ValueError(msg)
        if self.start < 0 or self.end < 0:
            msg = "primer binding coordinates must be non-negative"
            raise ValueError(msg)
        if self.mismatches < 0:
            msg = "primer binding mismatches must be non-negative"
            raise ValueError(msg)
        object.__setattr__(self, "primer_sequence", self.primer_sequence.upper())
        object.__setattr__(self, "binding_sequence", self.binding_sequence.upper())
        object.__setattr__(self, "tail_sequence", self.tail_sequence.upper())

    @property
    def binding_length(self) -> int:
        """Return the annealing length."""

        return len(self.binding_sequence)


@dataclass(frozen=True, slots=True)
class DesignedPrimerPair:
    """Primer pair returned by primer3."""

    forward_sequence: str
    reverse_sequence: str
    forward_start: int
    forward_length: int
    reverse_start: int
    reverse_length: int
    product_size: int
    forward_tm: float | None = None
    reverse_tm: float | None = None
    forward_gc_percent: float | None = None
    reverse_gc_percent: float | None = None
    penalty: float | None = None


@dataclass(frozen=True, slots=True)
class PCRProduct:
    """A simulated PCR product."""

    sequence: str
    start: int
    end: int
    template_length: int
    forward_binding: PrimerBinding
    reverse_binding: PrimerBinding
    wraps_origin: bool = False
    full_circle: bool = False
    provenance: ProvenanceEvent | None = None

    @property
    def length(self) -> int:
        """Return product sequence length, including primer tails."""

        return len(self.sequence)

    @property
    def coordinates(self) -> tuple[int, int]:
        """Return template start/end coordinates for the amplified interval."""

        return (self.start, self.end)


def find_primer_bindings(
    record: SequenceRecord,
    primer_sequence: str,
    *,
    mismatches: int = 0,
    include_forward: bool = True,
    include_reverse: bool = True,
) -> tuple[PrimerBinding, ...]:
    """Find full-length primer binding sites on a template."""

    _require_dna(record)
    primer = _normalize_primer(primer_sequence)
    if mismatches < 0:
        msg = "mismatches must be non-negative"
        raise ValueError(msg)

    bindings: list[PrimerBinding] = []
    if include_forward:
        bindings.extend(
            _find_matches(
                record,
                primer,
                primer,
                "",
                1,
                mismatches,
                annealing_tm=_tm_for_sequence(primer),
            )
        )
    if include_reverse:
        bindings.extend(
            _find_matches(
                record,
                _reverse_complement(primer),
                primer,
                "",
                -1,
                mismatches,
                annealing_tm=_tm_for_sequence(primer),
            )
        )
    annotated = _annotate_binding_uniqueness(bindings)
    return tuple(sorted(annotated, key=lambda binding: (binding.start, binding.strand, binding.mismatches)))


def primer_metrics(primer_sequence: str) -> PrimerMetrics:
    """Calculate basic primer metrics and primer3 metrics when available."""

    sequence = _normalize_primer(primer_sequence)
    gc_count = sequence.count("G") + sequence.count("C")
    at_count = sequence.count("A") + sequence.count("T")
    length = len(sequence)
    gc_percent = (gc_count / length * 100.0) if length else 0.0
    wallace_tm = float(2 * at_count + 4 * gc_count)
    primer3_tm: float | None = None
    hairpin = SecondaryStructureCheck(available=False)
    self_dimer = SecondaryStructureCheck(available=False)

    if primer3_bindings is not None:
        try:
            primer3_tm = float(primer3_bindings.calc_tm(sequence))
            hairpin = _secondary_structure_check(primer3_bindings.calc_hairpin(sequence))
            self_dimer = _secondary_structure_check(primer3_bindings.calc_homodimer(sequence))
        except Exception:
            hairpin = SecondaryStructureCheck(available=False)
            self_dimer = SecondaryStructureCheck(available=False)

    return PrimerMetrics(
        sequence=sequence,
        length=length,
        gc_percent=gc_percent,
        wallace_tm=wallace_tm,
        primer3_tm=primer3_tm,
        hairpin=hairpin,
        self_dimer=self_dimer,
    )


def design_primers(
    record: SequenceRecord,
    *,
    target_region: tuple[int, int],
    product_size_range: tuple[int, int] | Iterable[tuple[int, int]],
    primer_length_range: tuple[int, int] = (18, 25),
    tm_range: tuple[float, float] = (57.0, 63.0),
    gc_range: tuple[float, float] = (40.0, 60.0),
    num_return: int = 5,
) -> tuple[DesignedPrimerPair, ...]:
    """Design primer pairs using primer3-py."""

    _require_dna(record)
    if primer3_bindings is None:
        msg = "primer3-py is required for primer design"
        raise Primer3UnavailableError(msg)
    start, end = _validate_range(target_region, record.length, "target region")
    min_size, max_size = _validate_ordered_pair(product_size_range, "product size range")
    min_len, max_len = _validate_ordered_pair(primer_length_range, "primer length range")
    min_tm, max_tm = _validate_ordered_pair(tm_range, "Tm range")
    min_gc, max_gc = _validate_ordered_pair(gc_range, "GC range")
    opt_len = round((min_len + max_len) / 2)
    opt_tm = (min_tm + max_tm) / 2

    result = primer3_bindings.design_primers(
        {
            "SEQUENCE_ID": record.id,
            "SEQUENCE_TEMPLATE": record.sequence,
            "SEQUENCE_TARGET": [start, end - start],
        },
        {
            "PRIMER_PICK_LEFT_PRIMER": 1,
            "PRIMER_PICK_RIGHT_PRIMER": 1,
            "PRIMER_PRODUCT_SIZE_RANGE": [[min_size, max_size]],
            "PRIMER_MIN_SIZE": min_len,
            "PRIMER_OPT_SIZE": opt_len,
            "PRIMER_MAX_SIZE": max_len,
            "PRIMER_MIN_TM": min_tm,
            "PRIMER_OPT_TM": opt_tm,
            "PRIMER_MAX_TM": max_tm,
            "PRIMER_MIN_GC": min_gc,
            "PRIMER_MAX_GC": max_gc,
            "PRIMER_NUM_RETURN": num_return,
        },
    )
    return _designed_pairs_from_primer3(result)


def simulate_pcr(
    record: SequenceRecord,
    forward_primer: str,
    reverse_primer: str,
    *,
    mismatches: int = 0,
    min_anneal_length: int = 12,
    max_anneal_length: int | None = None,
    target_region: tuple[int, int] | None = None,
    pair_index: int | None = None,
    forward_binding_index: int | None = None,
    reverse_binding_index: int | None = None,
) -> PCRProduct:
    """Simulate PCR using a forward and reverse primer.

    Primer tails are inferred as non-template 5-prime prefixes. The 3-prime suffix of
    each primer must anneal to the template.
    """

    _require_dna(record)
    validated_target = _validate_pcr_target_region(target_region, record.length, record.topology)
    forward_bindings = _find_annealing_bindings(
        record,
        forward_primer,
        expected_strand=1,
        mismatches=mismatches,
        min_anneal_length=min_anneal_length,
        max_anneal_length=max_anneal_length,
        target_region=validated_target,
    )
    reverse_bindings = _find_annealing_bindings(
        record,
        reverse_primer,
        expected_strand=-1,
        mismatches=mismatches,
        min_anneal_length=min_anneal_length,
        max_anneal_length=max_anneal_length,
        target_region=validated_target,
    )
    if forward_binding_index is not None:
        forward_bindings = (_select_binding(forward_bindings, forward_binding_index, "forward"),)
    if reverse_binding_index is not None:
        reverse_bindings = (_select_binding(reverse_bindings, reverse_binding_index, "reverse"),)

    products = tuple(
        product
        for forward_binding in forward_bindings
        for reverse_binding in reverse_bindings
        if (
            product := _pcr_product_from_bindings(record, forward_binding, reverse_binding)
        )
        is not None
    )
    products = tuple(sorted(products, key=lambda product: _product_sort_key(product, validated_target)))
    if validated_target is not None and products:
        products = _select_products_for_target(products, validated_target)
    if not products:
        if validated_target is not None:
            msg = "primers do not produce a PCR product spanning the requested target_region"
        else:
            msg = "primers do not face each other or do not bind the template in PCR-compatible orientations"
        raise PCRSimulationError(msg)
    if pair_index is not None:
        if pair_index < 0 or pair_index >= len(products):
            msg = f"pair_index {pair_index} is outside the {len(products)} possible PCR products"
            raise PCRSimulationError(msg)
        return products[pair_index]
    if len(products) > 1:
        msg = "multiple possible PCR products; select a primer binding pair"
        raise PCRAmbiguousProductError(msg)
    return products[0]


def _find_annealing_bindings(
    record: SequenceRecord,
    primer_sequence: str,
    *,
    expected_strand: int,
    mismatches: int,
    min_anneal_length: int,
    max_anneal_length: int | None = None,
    target_region: tuple[int, int] | None = None,
) -> tuple[PrimerBinding, ...]:
    primer = _normalize_primer(primer_sequence)
    if min_anneal_length <= 0:
        msg = "min_anneal_length must be positive"
        raise ValueError(msg)
    minimum = min(min_anneal_length, len(primer))
    maximum = len(primer) if max_anneal_length is None else min(max_anneal_length, len(primer))
    if maximum < minimum:
        msg = "max_anneal_length must be greater than or equal to min_anneal_length"
        raise ValueError(msg)

    candidates: list[PrimerBinding] = []
    for anneal_length in range(maximum, minimum - 1, -1):
        binding_sequence = primer[-anneal_length:]
        tail_sequence = primer[: len(primer) - anneal_length]
        template_pattern = (
            binding_sequence if expected_strand == 1 else _reverse_complement(binding_sequence)
        )
        bindings = _find_matches(
            record,
            template_pattern,
            primer,
            tail_sequence,
            expected_strand,
            mismatches,
            annealing_tm=_tm_for_sequence(binding_sequence),
        )
        if bindings:
            candidates.extend(bindings)
    if not candidates:
        return ()

    deduplicated = _deduplicate_annealing_candidates(candidates, record.length, target_region)
    annotated = _annotate_binding_uniqueness(deduplicated)
    return tuple(sorted(annotated, key=lambda binding: _binding_sort_key(binding, target_region)))


def _find_matches(
    record: SequenceRecord,
    template_pattern: str,
    primer_sequence: str,
    tail_sequence: str,
    strand: int,
    mismatch_tolerance: int,
    *,
    annealing_tm: float | None = None,
) -> list[PrimerBinding]:
    sequence_length = record.length
    pattern_length = len(template_pattern)
    if sequence_length == 0 or pattern_length == 0:
        return []
    if record.topology is MoleculeTopology.LINEAR and pattern_length > sequence_length:
        return []

    max_start = sequence_length if record.topology is MoleculeTopology.CIRCULAR else sequence_length - pattern_length + 1
    matches: list[PrimerBinding] = []
    seen: set[tuple[int, int, int]] = set()
    for start in range(max_start):
        template_window = _template_window(record.sequence, start, pattern_length, record.topology)
        if template_window is None:
            continue
        mismatch_count = _mismatch_count(template_window, template_pattern)
        if mismatch_count > mismatch_tolerance:
            continue
        raw_end = start + pattern_length
        wraps_origin = record.topology is MoleculeTopology.CIRCULAR and raw_end > sequence_length
        end = raw_end % sequence_length if wraps_origin else raw_end
        key = (start, end, strand)
        if key in seen:
            continue
        matches.append(
            PrimerBinding(
                primer_sequence=primer_sequence,
                binding_sequence=primer_sequence[len(tail_sequence) :],
                tail_sequence=tail_sequence,
                start=start,
                end=end,
                strand=strand,
                mismatches=mismatch_count,
                wraps_origin=wraps_origin,
                annealing_tm=annealing_tm,
            )
        )
        seen.add(key)
    return matches


def _deduplicate_annealing_candidates(
    candidates: Iterable[PrimerBinding],
    template_length: int,
    target_region: tuple[int, int] | None,
) -> tuple[PrimerBinding, ...]:
    """Keep the best 3-prime-anchored candidate for each physical binding site."""

    best_by_anchor: dict[tuple[int, int], PrimerBinding] = {}
    for binding in candidates:
        key = _binding_anchor_key(binding, template_length)
        previous = best_by_anchor.get(key)
        if previous is None or _binding_preference_key(binding, target_region) < _binding_preference_key(
            previous,
            target_region,
        ):
            best_by_anchor[key] = binding
    return tuple(best_by_anchor.values())


def _annotate_binding_uniqueness(bindings: Iterable[PrimerBinding]) -> tuple[PrimerBinding, ...]:
    binding_tuple = tuple(bindings)
    counts = Counter((binding.binding_sequence, binding.strand) for binding in binding_tuple)
    return tuple(
        _binding_with_uniqueness(
            binding,
            is_unique=counts[(binding.binding_sequence, binding.strand)] == 1,
        )
        for binding in binding_tuple
    )


def _binding_with_uniqueness(binding: PrimerBinding, *, is_unique: bool) -> PrimerBinding:
    return PrimerBinding(
        primer_sequence=binding.primer_sequence,
        binding_sequence=binding.binding_sequence,
        tail_sequence=binding.tail_sequence,
        start=binding.start,
        end=binding.end,
        strand=binding.strand,
        mismatches=binding.mismatches,
        wraps_origin=binding.wraps_origin,
        annealing_tm=binding.annealing_tm,
        is_unique=is_unique,
    )


def _binding_anchor_key(binding: PrimerBinding, template_length: int) -> tuple[int, int]:
    if binding.strand == 1:
        return (binding.strand, binding.end % template_length)
    return (binding.strand, binding.start % template_length)


def _binding_preference_key(
    binding: PrimerBinding,
    target_region: tuple[int, int] | None,
) -> tuple[int, int, int, float, int, int, int]:
    target_miss = 0 if target_region is not None and _binding_overlaps_target(binding, target_region) else 1
    return (
        target_miss,
        -binding.binding_length,
        binding.mismatches,
        -(binding.annealing_tm or 0.0),
        0 if binding.is_unique else 1,
        binding.start,
        binding.end,
    )


def _binding_sort_key(
    binding: PrimerBinding,
    target_region: tuple[int, int] | None,
) -> tuple[int, int, int, float, int, int, int]:
    return _binding_preference_key(binding, target_region)


def _binding_overlaps_target(binding: PrimerBinding, target_region: tuple[int, int]) -> bool:
    target_start, target_end = target_region
    return binding.start < target_end and target_start < binding.end


def _pcr_product_from_bindings(
    record: SequenceRecord,
    forward_binding: PrimerBinding,
    reverse_binding: PrimerBinding,
) -> PCRProduct | None:
    if forward_binding.strand != 1 or reverse_binding.strand != -1:
        return None

    if record.topology is MoleculeTopology.LINEAR:
        if forward_binding.end > reverse_binding.start:
            return None
        template_sequence = record.sequence[forward_binding.start : reverse_binding.end]
        wraps_origin = False
    else:
        distance_to_forward_end = _circular_distance(
            forward_binding.start,
            forward_binding.end,
            record.length,
        )
        distance_to_reverse_start = _circular_distance(
            forward_binding.start,
            reverse_binding.start,
            record.length,
        )
        distance_to_reverse_end = _circular_distance(
            forward_binding.start,
            reverse_binding.end,
            record.length,
        )
        full_circle = (
            reverse_binding.end == forward_binding.start
            and distance_to_reverse_start > 0
        )
        if full_circle:
            distance_to_reverse_end = record.length
        if not (0 < distance_to_forward_end <= distance_to_reverse_start < distance_to_reverse_end):
            return None
        template_sequence = _circular_subsequence(
            record.sequence,
            forward_binding.start,
            reverse_binding.end,
            full_circle=full_circle,
        )
        wraps_origin = full_circle or forward_binding.start > reverse_binding.end

    product_sequence = (
        forward_binding.tail_sequence
        + template_sequence
        + _reverse_complement(reverse_binding.tail_sequence)
    )
    output_record_id = f"{record.id}_pcr_product"
    provenance = ProvenanceEvent(
        operation="pcr",
        input_record_id=record.id,
        input_record_ids=(record.id,),
        input_version_id=record.version_id,
        input_version_ids=(record.version_id,),
        output_record_id=output_record_id,
        affected_ranges=_pcr_affected_ranges(
            record,
            forward_binding.start,
            reverse_binding.end,
            full_circle=record.topology is MoleculeTopology.CIRCULAR
            and len(template_sequence) == record.length,
        ),
        description=f"PCR amplified {record.id}",
        parameters={
            "template_id": record.id,
            "start": forward_binding.start,
            "end": reverse_binding.end,
            "length": len(product_sequence),
            "wraps_origin": wraps_origin,
            "full_circle": record.topology is MoleculeTopology.CIRCULAR
            and len(template_sequence) == record.length,
            "template_interval_semantics": (
                "full_circle"
                if record.topology is MoleculeTopology.CIRCULAR
                and len(template_sequence) == record.length
                else "half_open"
            ),
            "forward_primer": forward_binding.primer_sequence,
            "reverse_primer": reverse_binding.primer_sequence,
            "forward_binding": _binding_report(forward_binding),
            "reverse_binding": _binding_report(reverse_binding),
        },
    )
    provenance = provenance.with_output(
        output_record_id,
        output_version_id=_derive_version_id(
            record_id=output_record_id,
            input_version_id=record.version_id,
            event=provenance,
            sequence=product_sequence,
            topology=MoleculeTopology.LINEAR,
        ),
    )
    return PCRProduct(
        sequence=product_sequence,
        start=forward_binding.start,
        end=reverse_binding.end,
        template_length=record.length,
        forward_binding=forward_binding,
        reverse_binding=reverse_binding,
        wraps_origin=wraps_origin,
        full_circle=record.topology is MoleculeTopology.CIRCULAR and len(template_sequence) == record.length,
        provenance=provenance,
    )


def _template_window(
    sequence: str,
    start: int,
    length: int,
    topology: MoleculeTopology,
) -> str | None:
    end = start + length
    if topology is MoleculeTopology.LINEAR:
        return sequence[start:end] if end <= len(sequence) else None
    if end <= len(sequence):
        return sequence[start:end]
    return sequence[start:] + sequence[: end % len(sequence)]


def _mismatch_count(template_window: str, template_pattern: str) -> int:
    return sum(
        0 if IUPAC_DNA[template_base] & IUPAC_DNA[pattern_base] else 1
        for template_base, pattern_base in zip(template_window, template_pattern, strict=True)
    )


def _tm_for_sequence(sequence: str) -> float:
    if primer3_bindings is not None:
        try:
            return float(primer3_bindings.calc_tm(sequence))
        except Exception:
            pass
    gc_count = sequence.count("G") + sequence.count("C")
    at_count = sequence.count("A") + sequence.count("T")
    return float(2 * at_count + 4 * gc_count)


def _normalize_primer(sequence: str) -> str:
    if not isinstance(sequence, str):
        msg = "primer sequence must be a string"
        raise TypeError(msg)
    normalized = sequence.upper().replace("U", "T")
    if not normalized:
        msg = "primer sequence must not be empty"
        raise ValueError(msg)
    invalid = sorted(set(normalized) - set(IUPAC_DNA))
    if invalid:
        msg = f"invalid primer sequence characters: {''.join(invalid)}"
        raise ValueError(msg)
    return normalized


def _reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(DNA_COMPLEMENT)[::-1]


def _secondary_structure_check(result: Any) -> SecondaryStructureCheck:
    return SecondaryStructureCheck(
        available=True,
        structure_found=bool(result.structure_found),
        tm=float(result.tm),
        dg=float(result.dg),
    )


def _designed_pairs_from_primer3(result: dict[str, Any]) -> tuple[DesignedPrimerPair, ...]:
    pair_count = int(result.get("PRIMER_PAIR_NUM_RETURNED", 0))
    pairs: list[DesignedPrimerPair] = []
    for index in range(pair_count):
        left = result.get("PRIMER_LEFT", [{}] * pair_count)[index]
        right = result.get("PRIMER_RIGHT", [{}] * pair_count)[index]
        pair = result.get("PRIMER_PAIR", [{}] * pair_count)[index]
        left_coords = left.get("COORDS", result.get(f"PRIMER_LEFT_{index}"))
        right_coords = right.get("COORDS", result.get(f"PRIMER_RIGHT_{index}"))
        if not left_coords or not right_coords:
            continue
        pairs.append(
            DesignedPrimerPair(
                forward_sequence=_primer3_value(result, left, "SEQUENCE", f"PRIMER_LEFT_{index}_SEQUENCE"),
                reverse_sequence=_primer3_value(result, right, "SEQUENCE", f"PRIMER_RIGHT_{index}_SEQUENCE"),
                forward_start=int(left_coords[0]),
                forward_length=int(left_coords[1]),
                reverse_start=int(right_coords[0]) - int(right_coords[1]) + 1,
                reverse_length=int(right_coords[1]),
                product_size=int(
                    _primer3_value(result, pair, "PRODUCT_SIZE", f"PRIMER_PAIR_{index}_PRODUCT_SIZE")
                ),
                forward_tm=_optional_float(left.get("TM", result.get(f"PRIMER_LEFT_{index}_TM"))),
                reverse_tm=_optional_float(right.get("TM", result.get(f"PRIMER_RIGHT_{index}_TM"))),
                forward_gc_percent=_optional_float(
                    left.get("GC_PERCENT", result.get(f"PRIMER_LEFT_{index}_GC_PERCENT"))
                ),
                reverse_gc_percent=_optional_float(
                    right.get("GC_PERCENT", result.get(f"PRIMER_RIGHT_{index}_GC_PERCENT"))
                ),
                penalty=_optional_float(pair.get("PENALTY", result.get(f"PRIMER_PAIR_{index}_PENALTY"))),
            )
        )
    return tuple(pairs)


def _optional_float(value: Any) -> float | None:
    return None if value is None else float(value)


def _primer3_value(result: dict[str, Any], nested: dict[str, Any], nested_key: str, flat_key: str) -> Any:
    if nested_key in nested:
        return nested[nested_key]
    return result[flat_key]


def _validate_range(value: tuple[int, int], length: int, field_name: str) -> tuple[int, int]:
    start, end = value
    if start < 0 or end < start or end > length:
        msg = f"{field_name} must be a zero-based half-open interval within the template"
        raise ValueError(msg)
    return start, end


def _validate_ordered_pair(value: Any, field_name: str) -> tuple[Any, Any]:
    if isinstance(value, tuple) and len(value) == 2 and not isinstance(value[0], tuple):
        lower, upper = value
    else:
        ranges = tuple(value)
        if len(ranges) != 1:
            msg = f"{field_name} currently accepts one range"
            raise ValueError(msg)
        lower, upper = ranges[0]
    if lower > upper:
        msg = f"{field_name} minimum must be less than or equal to maximum"
        raise ValueError(msg)
    return lower, upper


def _select_binding(
    bindings: tuple[PrimerBinding, ...],
    index: int,
    label: str,
) -> PrimerBinding:
    if index < 0 or index >= len(bindings):
        msg = f"{label} binding index {index} is outside the {len(bindings)} available bindings"
        raise PCRSimulationError(msg)
    return bindings[index]


def _validate_pcr_target_region(
    target_region: tuple[int, int] | None,
    template_length: int,
    topology: MoleculeTopology,
) -> tuple[int, int] | None:
    if target_region is None:
        return None
    start, end = target_region
    if start < 0 or end < 0 or start > template_length or end > template_length:
        msg = "target_region must be within the template"
        raise ValueError(msg)
    if start >= end:
        if topology is MoleculeTopology.CIRCULAR:
            msg = "target_region must be a non-wrapping zero-based half-open interval"
        else:
            msg = "target_region start must be less than end"
        raise ValueError(msg)
    return (start, end)


def _product_sort_key(
    product: PCRProduct,
    target_region: tuple[int, int] | None,
) -> tuple[int, int, int, int, float, int, int, int]:
    target_excess = _product_target_excess(product, target_region) if target_region else None
    return (
        1 if target_region is not None and target_excess is None else 0,
        target_excess if target_excess is not None else 0,
        -(product.forward_binding.binding_length + product.reverse_binding.binding_length),
        product.forward_binding.mismatches + product.reverse_binding.mismatches,
        -((product.forward_binding.annealing_tm or 0.0) + (product.reverse_binding.annealing_tm or 0.0)),
        0 if product.forward_binding.is_unique and product.reverse_binding.is_unique else 1,
        product.start,
        product.end,
    )


def _select_products_for_target(
    products: tuple[PCRProduct, ...],
    target_region: tuple[int, int],
) -> tuple[PCRProduct, ...]:
    covering = tuple(product for product in products if _product_target_excess(product, target_region) is not None)
    if not covering:
        return ()
    best_excess = min(_product_target_excess(product, target_region) for product in covering)
    return tuple(
        product
        for product in covering
        if _product_target_excess(product, target_region) == best_excess
    )


def _product_target_excess(product: PCRProduct, target_region: tuple[int, int] | None) -> int | None:
    if target_region is None:
        return 0
    if not _product_covers_target(product, target_region):
        return None
    return _template_span_length(product) - (target_region[1] - target_region[0])


def _product_covers_target(product: PCRProduct, target_region: tuple[int, int]) -> bool:
    target_start, target_end = target_region
    if product.full_circle:
        return True
    for product_start, product_end in _product_intervals(product):
        if product_start <= target_start and target_end <= product_end:
            return True
    return False


def _product_intervals(product: PCRProduct) -> tuple[tuple[int, int], ...]:
    if product.full_circle:
        return ((0, product.template_length),)
    if not product.wraps_origin:
        return ((product.start, product.end),)
    return ((product.start, product.template_length), (0, product.end))


def _template_span_length(product: PCRProduct) -> int:
    if product.full_circle:
        return product.template_length
    if not product.wraps_origin:
        return product.end - product.start
    return product.template_length - product.start + product.end


def _binding_report(binding: PrimerBinding) -> dict[str, Any]:
    return {
        "start": binding.start,
        "end": binding.end,
        "strand": binding.strand,
        "anneal_length": binding.binding_length,
        "tail_length": len(binding.tail_sequence),
        "mismatches": binding.mismatches,
        "annealing_tm": binding.annealing_tm,
        "unique": binding.is_unique,
        "wraps_origin": binding.wraps_origin,
    }


def _circular_distance(start: int, end: int, length: int) -> int:
    return end - start if end >= start else length - start + end


def _circular_subsequence(
    sequence: str,
    start: int,
    end: int,
    *,
    full_circle: bool = False,
) -> str:
    if full_circle:
        return sequence[start:] + sequence[:start]
    if start < end:
        return sequence[start:end]
    return sequence[start:] + sequence[:end]


def _pcr_affected_ranges(
    record: SequenceRecord,
    start: int,
    end: int,
    *,
    full_circle: bool = False,
) -> tuple[tuple[int, int], ...]:
    if record.topology is MoleculeTopology.CIRCULAR and full_circle:
        if start == 0:
            return ((0, record.length),)
        return ((start, record.length), (0, start))
    if record.topology is MoleculeTopology.LINEAR or start <= end:
        return ((start, end),)
    return ((start, record.length), (0, end))


def _require_dna(record: SequenceRecord) -> None:
    if record.molecule_type is not MoleculeType.DNA:
        msg = "primer analysis requires a DNA sequence"
        raise ValueError(msg)
