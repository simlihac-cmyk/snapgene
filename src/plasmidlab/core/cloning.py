"""Core cloning simulation algorithms."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Any

from plasmidlab.core.models import (
    Feature,
    FeatureSegment,
    MoleculeTopology,
    MoleculeType,
    ProvenanceEvent,
    SequenceRecord,
    _derive_version_id,
)
from plasmidlab.core.restriction import (
    FragmentEnd,
    FragmentEndSide,
    RestrictionEnzyme,
    RestrictionSite,
    compatible_fragment_ends,
    find_restriction_sites,
    fragment_end_for_site,
    resolve_enzymes,
)


DNA_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


class CloningError(ValueError):
    """Raised when a cloning simulation cannot produce a product."""


@dataclass(frozen=True, slots=True)
class AssemblyOverlap:
    """Overlap used between two assembly fragments."""

    left_fragment_index: int
    right_fragment_index: int
    left_start: int
    left_end: int
    right_start: int
    right_end: int
    length: int
    sequence: str


@dataclass(frozen=True, slots=True)
class CloningResult:
    """Result from a cloning or assembly operation."""

    product: SequenceRecord
    warnings: tuple[str, ...] = ()
    overlaps: tuple[AssemblyOverlap, ...] = ()

    def __getattr__(self, name: str) -> Any:
        """Delegate unknown attributes to the product record for convenience."""

        return getattr(self.product, name)


@dataclass(frozen=True, slots=True)
class _Fragment:
    sequence: str
    features: tuple[Feature, ...]
    left_end: FragmentEnd
    right_end: FragmentEnd


@dataclass(frozen=True, slots=True)
class _GoldenGateFragment:
    core_sequence: str
    features: tuple[Feature, ...]
    left_end: FragmentEnd
    right_end: FragmentEnd

    @property
    def left_overhang(self) -> str:
        return self.left_end.top_strand_overhang_sequence

    @property
    def right_overhang(self) -> str:
        return self.right_end.top_strand_overhang_sequence


def restriction_clone(
    vector: SequenceRecord,
    insert: SequenceRecord,
    vector_enzymes: Any | Iterable[Any],
    insert_enzymes: Any | Iterable[Any] | None = None,
    *,
    insert_orientation: str = "forward",
    allow_ambiguous_overhangs: bool = False,
    product_id: str | None = None,
) -> CloningResult:
    """Cut vector and insert with restriction enzymes and ligate compatible ends."""

    _require_dna(vector)
    _require_dna(insert)
    vector_enzyme_set = _one_or_two_enzymes(vector_enzymes, "vector enzymes")
    insert_enzyme_set = (
        vector_enzyme_set if insert_enzymes is None else _one_or_two_enzymes(insert_enzymes, "insert enzymes")
    )
    vector_fragment = _restriction_vector_fragment(vector, vector_enzyme_set)
    insert_fragment = _restriction_insert_fragment(insert, insert_enzyme_set)
    if insert_orientation == "reverse":
        insert_fragment = _reverse_fragment(insert_fragment)
    elif insert_orientation != "forward":
        msg = "insert_orientation must be 'forward' or 'reverse'"
        raise ValueError(msg)

    _require_compatible(
        vector_fragment.right_end,
        insert_fragment.left_end,
        "vector right",
        "insert left",
        allow_ambiguous=allow_ambiguous_overhangs,
    )
    _require_compatible(
        insert_fragment.right_end,
        vector_fragment.left_end,
        "insert right",
        "vector left",
        allow_ambiguous=allow_ambiguous_overhangs,
    )

    product_sequence = vector_fragment.sequence + insert_fragment.sequence
    features = vector_fragment.features + _shift_features(insert_fragment.features, len(vector_fragment.sequence))
    resolved_product_id = product_id or f"{vector.id}_{insert.id}_restriction_clone"
    product = _new_product_record(
        id=resolved_product_id,
        sequence=product_sequence,
        topology=MoleculeTopology.CIRCULAR,
        features=features,
        histories=(vector.history, insert.history),
        event=ProvenanceEvent(
            operation="restriction_clone",
            input_record_id=vector.id,
            input_record_ids=(vector.id, insert.id),
            input_version_id=vector.version_id,
            input_version_ids=(vector.version_id, insert.version_id),
            output_record_id=resolved_product_id,
            affected_ranges=((0, len(product_sequence)),),
            description=f"Restriction cloned {insert.id} into {vector.id}",
            parameters={
                "vector_id": vector.id,
                "insert_id": insert.id,
                "vector_enzymes": tuple(enzyme.name for enzyme in vector_enzyme_set),
                "insert_enzymes": tuple(enzyme.name for enzyme in insert_enzyme_set),
                "insert_orientation": insert_orientation,
                "allow_ambiguous_overhangs": allow_ambiguous_overhangs,
            },
        ),
    )
    return CloningResult(product=product)


def gibson_assembly(
    fragments: Sequence[SequenceRecord],
    *,
    min_overlap: int = 20,
    circular: bool = False,
    product_id: str = "gibson_assembly",
) -> CloningResult:
    """Assemble ordered fragments by terminal sequence overlaps."""

    if len(fragments) < 2:
        msg = "Gibson assembly requires at least two fragments"
        raise CloningError(msg)
    if min_overlap <= 0:
        msg = "minimum overlap length must be positive"
        raise ValueError(msg)
    for fragment in fragments:
        _require_dna(fragment)

    overlaps: list[AssemblyOverlap] = []
    internal_overlap_lengths: list[int] = []
    warnings: list[str] = []
    for index, (left, right) in enumerate(zip(fragments, fragments[1:], strict=False)):
        overlap = _terminal_overlap(left.sequence, right.sequence, min_overlap)
        if overlap == 0:
            msg = f"fragments {index} and {index + 1} do not share a terminal overlap"
            raise CloningError(msg)
        internal_overlap_lengths.append(overlap)
        overlaps.append(
            AssemblyOverlap(
                left_fragment_index=index,
                right_fragment_index=index + 1,
                left_start=len(left.sequence) - overlap,
                left_end=len(left.sequence),
                right_start=0,
                right_end=overlap,
                length=overlap,
                sequence=left.sequence[-overlap:],
            )
        )

    circular_overlap = _terminal_overlap(fragments[-1].sequence, fragments[0].sequence, min_overlap)
    if circular:
        if circular_overlap == 0:
            msg = "last and first fragments do not share the requested circular overlap"
            raise CloningError(msg)
        overlaps.append(
            AssemblyOverlap(
                left_fragment_index=len(fragments) - 1,
                right_fragment_index=0,
                left_start=len(fragments[-1].sequence) - circular_overlap,
                left_end=len(fragments[-1].sequence),
                right_start=0,
                right_end=circular_overlap,
                length=circular_overlap,
                sequence=fragments[-1].sequence[-circular_overlap:],
            )
        )
    elif circular_overlap:
        warnings.append("terminal fragments also share an overlap; circular assembly may be possible")

    sequence_parts = [fragments[0].sequence]
    feature_parts = [_project_features(fragments[0], [(0, fragments[0].length)], 0)]
    offset = len(sequence_parts[0])
    for index, fragment in enumerate(fragments[1:], start=1):
        overlap = internal_overlap_lengths[index - 1]
        added_sequence = fragment.sequence[overlap:]
        sequence_parts.append(added_sequence)
        feature_parts.append(_project_features(fragment, [(overlap, fragment.length)], offset))
        offset += len(added_sequence)

    assembled_sequence = "".join(sequence_parts)
    features = tuple(feature for part in feature_parts for feature in part)
    if circular:
        assembled_sequence = assembled_sequence[:-circular_overlap]
        features = _trim_features_at_end(features, len(assembled_sequence))

    product = _new_product_record(
        id=product_id,
        sequence=assembled_sequence,
        topology=MoleculeTopology.CIRCULAR if circular else MoleculeTopology.LINEAR,
        features=features,
        histories=tuple(fragment.history for fragment in fragments),
        event=ProvenanceEvent(
            operation="gibson_assembly",
            input_record_id=fragments[0].id,
            input_record_ids=tuple(fragment.id for fragment in fragments),
            input_version_id=fragments[0].version_id,
            input_version_ids=tuple(fragment.version_id for fragment in fragments),
            output_record_id=product_id,
            affected_ranges=((0, len(assembled_sequence)),),
            description=f"Gibson assembled {len(fragments)} fragments",
            parameters={
                "fragment_ids": tuple(fragment.id for fragment in fragments),
                "min_overlap": min_overlap,
                "circular": circular,
            },
        ),
    )
    return CloningResult(product=product, warnings=tuple(warnings), overlaps=tuple(overlaps))


def golden_gate_assembly(
    fragments: Sequence[SequenceRecord],
    enzyme: Any,
    *,
    circular: bool = True,
    allow_ambiguous_overhangs: bool = False,
    product_id: str = "golden_gate_assembly",
) -> CloningResult:
    """Assemble ordered Type IIS fragments by compatible generated overhangs."""

    if len(fragments) < 2:
        msg = "Golden Gate assembly requires at least two fragments"
        raise CloningError(msg)
    enzyme_definition = resolve_enzymes((enzyme,))[0]
    if enzyme_definition.top_cut_offset <= enzyme_definition.recognition_length:
        msg = "Golden Gate assembly requires a Type IIS enzyme cutting outside its recognition site"
        raise CloningError(msg)
    for fragment in fragments:
        _require_dna(fragment)

    prepared = tuple(_golden_gate_fragment(fragment, enzyme_definition) for fragment in fragments)
    warnings = list(_golden_gate_warnings(prepared))
    for index, (left, right) in enumerate(zip(prepared, prepared[1:], strict=False)):
        _require_compatible(
            left.right_end,
            right.left_end,
            f"fragment {index} right",
            f"fragment {index + 1} left",
            allow_ambiguous=allow_ambiguous_overhangs,
        )
    if circular:
        _require_compatible(
            prepared[-1].right_end,
            prepared[0].left_end,
            "last fragment right",
            "first fragment left",
            allow_ambiguous=allow_ambiguous_overhangs,
        )

    overlap_length = len(prepared[0].left_overhang)
    assembled = prepared[0].core_sequence
    features = list(prepared[0].features)
    offset = len(assembled)
    overlaps: list[AssemblyOverlap] = []
    for index, fragment in enumerate(prepared[1:], start=1):
        assembled += fragment.core_sequence[overlap_length:]
        features.extend(_shift_features(fragment.features, offset - overlap_length))
        overlaps.append(
            AssemblyOverlap(
                left_fragment_index=index - 1,
                right_fragment_index=index,
                left_start=len(prepared[index - 1].core_sequence) - overlap_length,
                left_end=len(prepared[index - 1].core_sequence),
                right_start=0,
                right_end=overlap_length,
                length=overlap_length,
                sequence=fragment.left_overhang,
            )
        )
        offset = len(assembled)

    if circular:
        assembled = assembled[:-overlap_length]
        features = list(_trim_features_at_end(tuple(features), len(assembled)))
        overlaps.append(
            AssemblyOverlap(
                left_fragment_index=len(prepared) - 1,
                right_fragment_index=0,
                left_start=len(prepared[-1].core_sequence) - overlap_length,
                left_end=len(prepared[-1].core_sequence),
                right_start=0,
                right_end=overlap_length,
                length=overlap_length,
                sequence=prepared[0].left_overhang,
            )
        )

    product = _new_product_record(
        id=product_id,
        sequence=assembled,
        topology=MoleculeTopology.CIRCULAR if circular else MoleculeTopology.LINEAR,
        features=tuple(features),
        histories=tuple(fragment.history for fragment in fragments),
        event=ProvenanceEvent(
            operation="golden_gate_assembly",
            input_record_id=fragments[0].id,
            input_record_ids=tuple(fragment.id for fragment in fragments),
            input_version_id=fragments[0].version_id,
            input_version_ids=tuple(fragment.version_id for fragment in fragments),
            output_record_id=product_id,
            affected_ranges=((0, len(assembled)),),
            description=f"Golden Gate assembled {len(fragments)} fragments with {enzyme_definition.name}",
            parameters={
                "fragment_ids": tuple(fragment.id for fragment in fragments),
                "enzyme": enzyme_definition.name,
                "circular": circular,
                "allow_ambiguous_overhangs": allow_ambiguous_overhangs,
            },
        ),
    )
    return CloningResult(product=product, warnings=tuple(warnings), overlaps=tuple(overlaps))


def inverse_pcr_mutagenesis(
    template: SequenceRecord,
    start: int,
    end: int,
    *,
    insertion_sequence: str = "",
    product_id: str | None = None,
) -> CloningResult:
    """Simulate inverse-PCR deletion/insertion mutagenesis and circularization."""

    _require_dna(template)
    edited = template.replace(start, end, insertion_sequence)
    resolved_product_id = product_id or f"{template.id}_inverse_pcr"
    product = _new_product_record(
        id=resolved_product_id,
        sequence=edited.sequence,
        topology=MoleculeTopology.CIRCULAR,
        features=edited.features,
        histories=(template.history,),
        event=ProvenanceEvent(
            operation="inverse_pcr_mutagenesis",
            input_record_id=template.id,
            input_record_ids=(template.id,),
            input_version_id=template.version_id,
            input_version_ids=(template.version_id,),
            output_record_id=resolved_product_id,
            affected_ranges=((start, end),),
            description=f"Inverse-PCR mutagenesis of {template.id}",
            parameters={
                "template_id": template.id,
                "start": start,
                "end": end,
                "inserted_length": len(insertion_sequence),
            },
        ),
    )
    return CloningResult(product=product)


def _restriction_vector_fragment(
    vector: SequenceRecord,
    enzymes: tuple[RestrictionEnzyme, ...],
) -> _Fragment:
    if vector.topology is not MoleculeTopology.CIRCULAR:
        msg = "restriction cloning currently expects a circular vector"
        raise CloningError(msg)
    if len(enzymes) == 1:
        site = _single_site(vector, enzymes[0])
        cut = site.top_cut
        return _Fragment(
            sequence=vector.sequence[cut:] + vector.sequence[:cut],
            features=_project_features(vector, [(cut, vector.length), (0, cut)], 0),
            left_end=fragment_end_for_site(vector, site, FragmentEndSide.LEFT),
            right_end=fragment_end_for_site(vector, site, FragmentEndSide.RIGHT),
        )

    first_site = _single_site(vector, enzymes[0])
    second_site = _single_site(vector, enzymes[1])
    intervals = _circular_intervals(second_site.top_cut, first_site.top_cut, vector.length)
    return _Fragment(
        sequence=_sequence_from_intervals(vector.sequence, intervals),
        features=_project_features(vector, intervals, 0),
        left_end=fragment_end_for_site(vector, second_site, FragmentEndSide.LEFT),
        right_end=fragment_end_for_site(vector, first_site, FragmentEndSide.RIGHT),
    )


def _restriction_insert_fragment(
    insert: SequenceRecord,
    enzymes: tuple[RestrictionEnzyme, ...],
) -> _Fragment:
    if len(enzymes) == 1:
        sites = _sites_for_enzyme(insert, enzymes[0])
        if len(sites) != 2:
            msg = "single-enzyme insert cloning requires exactly two insert cut sites"
            raise CloningError(msg)
        left_site, right_site = sorted(sites, key=lambda site: site.top_cut)
        intervals = [(left_site.top_cut, right_site.top_cut)]
        return _Fragment(
            sequence=_sequence_from_intervals(insert.sequence, intervals),
            features=_project_features(insert, intervals, 0),
            left_end=fragment_end_for_site(insert, left_site, FragmentEndSide.LEFT),
            right_end=fragment_end_for_site(insert, right_site, FragmentEndSide.RIGHT),
        )

    left_site = _single_site(insert, enzymes[0])
    right_site = _single_site(insert, enzymes[1])
    if left_site.top_cut >= right_site.top_cut:
        msg = "insert enzyme sites must be ordered on the linear insert"
        raise CloningError(msg)
    intervals = [(left_site.top_cut, right_site.top_cut)]
    return _Fragment(
        sequence=_sequence_from_intervals(insert.sequence, intervals),
        features=_project_features(insert, intervals, 0),
        left_end=fragment_end_for_site(insert, left_site, FragmentEndSide.LEFT),
        right_end=fragment_end_for_site(insert, right_site, FragmentEndSide.RIGHT),
    )


def _single_site(record: SequenceRecord, enzyme: RestrictionEnzyme) -> RestrictionSite:
    sites = _sites_for_enzyme(record, enzyme)
    if len(sites) != 1:
        msg = f"{enzyme.name} must cut {record.id} exactly once; found {len(sites)} sites"
        raise CloningError(msg)
    return sites[0]


def _sites_for_enzyme(record: SequenceRecord, enzyme: RestrictionEnzyme) -> tuple[RestrictionSite, ...]:
    return tuple(site for site in find_restriction_sites(record, (enzyme,)) if site.enzyme_name == enzyme.name)


def _golden_gate_fragment(
    record: SequenceRecord,
    enzyme: RestrictionEnzyme,
) -> _GoldenGateFragment:
    sites = find_restriction_sites(record, (enzyme,))
    plus_sites = [site for site in sites if site.strand == 1]
    minus_sites = [site for site in sites if site.strand == -1]
    if len(plus_sites) != 1 or len(minus_sites) != 1:
        msg = f"{record.id} must contain one forward and one reverse {enzyme.name} site"
        raise CloningError(msg)
    left_site = plus_sites[0]
    right_site = minus_sites[0]
    if left_site.top_cut >= right_site.top_cut:
        msg = f"{record.id} Type IIS sites do not face outward around an insert"
        raise CloningError(msg)
    left_start, left_end = sorted((left_site.top_cut, left_site.bottom_cut))
    right_start, right_end = sorted((right_site.top_cut, right_site.bottom_cut))
    if left_end > right_start:
        msg = f"{record.id} Type IIS cut positions overlap"
        raise CloningError(msg)

    return _GoldenGateFragment(
        core_sequence=record.sequence[left_start:right_end],
        features=_project_features(record, [(left_start, right_end)], 0),
        left_end=fragment_end_for_site(record, left_site, FragmentEndSide.LEFT),
        right_end=fragment_end_for_site(record, right_site, FragmentEndSide.RIGHT),
    )


def _golden_gate_warnings(fragments: tuple[_GoldenGateFragment, ...]) -> tuple[str, ...]:
    warnings: list[str] = []
    left_counts = Counter(fragment.left_overhang for fragment in fragments)
    right_counts = Counter(fragment.right_overhang for fragment in fragments)
    all_counts = Counter(
        overhang
        for fragment in fragments
        for overhang in (fragment.left_overhang, fragment.right_overhang)
    )
    duplicated_left = sorted(overhang for overhang, count in left_counts.items() if count > 1)
    duplicated_right = sorted(overhang for overhang, count in right_counts.items() if count > 1)
    ambiguous = sorted(
        overhang
        for overhang, count in all_counts.items()
        if count > 2 or left_counts[overhang] > 1 or right_counts[overhang] > 1
    )
    if duplicated_left:
        warnings.append(f"duplicated left overhangs: {', '.join(duplicated_left)}")
    if duplicated_right:
        warnings.append(f"duplicated right overhangs: {', '.join(duplicated_right)}")
    if ambiguous:
        warnings.append(f"ambiguous Golden Gate overhang reuse: {', '.join(ambiguous)}")
    return tuple(warnings)


def _terminal_overlap(left: str, right: str, min_overlap: int) -> int:
    max_overlap = min(len(left), len(right))
    for overlap in range(max_overlap, min_overlap - 1, -1):
        if left[-overlap:] == right[:overlap]:
            return overlap
    return 0


def _require_compatible(
    left: FragmentEnd,
    right: FragmentEnd,
    left_label: str,
    right_label: str,
    *,
    allow_ambiguous: bool = False,
) -> None:
    if compatible_fragment_ends(left, right, allow_ambiguous=allow_ambiguous):
        return
    msg = f"incompatible ends: {left_label} cannot ligate to {right_label}; overhang does not match"
    raise CloningError(msg)


def _reverse_fragment(fragment: _Fragment) -> _Fragment:
    sequence = _reverse_complement(fragment.sequence)
    return _Fragment(
        sequence=sequence,
        features=_reverse_features(fragment.features, len(fragment.sequence)),
        left_end=_reverse_fragment_end(fragment.right_end, FragmentEndSide.LEFT),
        right_end=_reverse_fragment_end(fragment.left_end, FragmentEndSide.RIGHT),
    )


def _reverse_fragment_end(end: FragmentEnd, side: FragmentEndSide) -> FragmentEnd:
    overhang_sequence = _reverse_complement(end.overhang_sequence) if end.overhang_sequence else ""
    top_strand_sequence = (
        _reverse_complement(end.top_strand_overhang_sequence)
        if end.top_strand_overhang_sequence
        else ""
    )
    return FragmentEnd(
        side=side,
        kind=end.kind,
        overhang_sequence=overhang_sequence,
        top_cut=end.top_cut,
        bottom_cut=end.bottom_cut,
        source_enzyme=end.source_enzyme,
        source_site_coordinate=end.source_site_coordinate,
        orientation=-end.orientation,
        top_strand_overhang_sequence=top_strand_sequence,
        cut_geometry=end.cut_geometry,
    )


def _reverse_features(features: Iterable[Feature], sequence_length: int) -> tuple[Feature, ...]:
    reversed_features: list[Feature] = []
    for feature in features:
        reversed_features.append(
            feature.with_segments(
                segment.reverse_complement(sequence_length) for segment in reversed(feature.segments)
            )
        )
    return tuple(reversed_features)


def _project_features(
    record: SequenceRecord,
    intervals: Iterable[tuple[int, int]],
    offset: int,
) -> tuple[Feature, ...]:
    interval_tuple = tuple(intervals)
    projected_features: list[Feature] = []
    for feature in record.features:
        segments: list[FeatureSegment] = []
        cursor = offset
        for interval_start, interval_end in interval_tuple:
            for segment in feature.segments:
                overlap_start = max(segment.start, interval_start)
                overlap_end = min(segment.end, interval_end)
                if overlap_start < overlap_end:
                    segments.append(
                        FeatureSegment(
                            cursor + overlap_start - interval_start,
                            cursor + overlap_end - interval_start,
                            segment.strand,
                        )
                    )
            cursor += interval_end - interval_start
        if segments:
            projected_features.append(feature.with_segments(_merge_segments(segments)))
    return tuple(projected_features)


def _shift_features(features: Iterable[Feature], offset: int) -> tuple[Feature, ...]:
    return tuple(
        feature.with_segments(segment.shift(offset) for segment in feature.segments)
        for feature in features
    )


def _trim_features_at_end(features: tuple[Feature, ...], new_length: int) -> tuple[Feature, ...]:
    trimmed: list[Feature] = []
    for feature in features:
        segments = [
            FeatureSegment(segment.start, min(segment.end, new_length), segment.strand)
            for segment in feature.segments
            if segment.start < new_length
        ]
        if segments:
            trimmed.append(feature.with_segments(_merge_segments(segments)))
    return tuple(trimmed)


def _merge_segments(segments: Iterable[FeatureSegment]) -> tuple[FeatureSegment, ...]:
    sorted_segments = sorted(segments, key=lambda segment: (segment.start, segment.end, segment.strand))
    if not sorted_segments:
        return ()
    merged = [sorted_segments[0]]
    for segment in sorted_segments[1:]:
        previous = merged[-1]
        if previous.strand == segment.strand and previous.end >= segment.start:
            merged[-1] = FeatureSegment(previous.start, max(previous.end, segment.end), previous.strand)
        else:
            merged.append(segment)
    return tuple(merged)


def _new_product_record(
    *,
    id: str,
    sequence: str,
    topology: MoleculeTopology,
    features: tuple[Feature, ...],
    histories: Iterable[tuple[ProvenanceEvent, ...]],
    event: ProvenanceEvent,
) -> SequenceRecord:
    output_version_id = event.output_version_id or _derive_version_id(
        record_id=id,
        input_version_id=event.input_version_id or id,
        event=event,
        sequence=sequence,
        topology=topology,
    )
    if event.output_record_id is None or event.output_version_id is None:
        event = event.with_output(id, output_version_id=output_version_id)
    history = tuple(event for history in histories for event in history) + (event,)
    return SequenceRecord(
        id=id,
        sequence=sequence,
        molecule_type=MoleculeType.DNA,
        topology=topology,
        version_id=output_version_id,
        features=features,
        history=history,
    )


def _one_or_two_enzymes(enzymes: Any | Iterable[Any], label: str) -> tuple[RestrictionEnzyme, ...]:
    resolved = resolve_enzymes(enzymes)
    if len(resolved) not in (1, 2):
        msg = f"{label} must contain one or two enzymes"
        raise ValueError(msg)
    return resolved


def _circular_intervals(start: int, end: int, length: int) -> tuple[tuple[int, int], ...]:
    if start == end:
        return ((0, length),)
    if start < end:
        return ((start, end),)
    return ((start, length), (0, end))


def _sequence_from_intervals(sequence: str, intervals: Iterable[tuple[int, int]]) -> str:
    return "".join(sequence[start:end] for start, end in intervals)


def _reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(DNA_COMPLEMENT)[::-1]


def _require_dna(record: SequenceRecord) -> None:
    if record.molecule_type is not MoleculeType.DNA:
        msg = "cloning simulations require DNA sequence records"
        raise ValueError(msg)
