"""Foundational biological data models and edit operations for PlasmidLab."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from datetime import UTC, datetime
from enum import StrEnum
from types import MappingProxyType
from typing import Any

QualifierValue = str | tuple[str, ...]


class MoleculeTopology(StrEnum):
    """Physical topology of a sequence molecule."""

    LINEAR = "linear"
    CIRCULAR = "circular"


# Backward-compatible short name used by the initial scaffold.
Topology = MoleculeTopology


class MoleculeType(StrEnum):
    """Supported biological molecule classes."""

    DNA = "DNA"
    RNA = "RNA"
    PROTEIN = "protein"


DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN")
RNA_ALPHABET = frozenset("ACGURYSWKMBDHVN")
PROTEIN_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWYBXZJUO*")

DNA_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
RNA_COMPLEMENT = str.maketrans("ACGURYSWKMBDHVN", "UGCAYRSWMKVHDBN")


def _validate_text(value: str, field_name: str) -> str:
    if not isinstance(value, str):
        msg = f"{field_name} must be a string"
        raise TypeError(msg)
    if not value:
        msg = f"{field_name} must not be empty"
        raise ValueError(msg)
    if value != value.strip():
        msg = f"{field_name} must not have leading or trailing whitespace"
        raise ValueError(msg)
    if any(ord(character) < 32 for character in value):
        msg = f"{field_name} must not contain control characters"
        raise ValueError(msg)
    return value


def _validate_optional_text(value: str | None, field_name: str) -> str | None:
    if value is None:
        return None
    return _validate_text(value, field_name)


def _validate_strand(strand: int) -> int:
    if strand not in (-1, 0, 1):
        msg = "strand must be +1, -1, or 0"
        raise ValueError(msg)
    return strand


def _freeze_qualifier_mapping(values: Mapping[str, Any] | None) -> Mapping[str, QualifierValue]:
    frozen: dict[str, QualifierValue] = {}
    for key, value in (values or {}).items():
        qualifier_key = _validate_text(str(key), "qualifier key")
        if isinstance(value, str):
            frozen[qualifier_key] = value
        elif isinstance(value, tuple | list):
            frozen[qualifier_key] = tuple(str(item) for item in value)
        else:
            frozen[qualifier_key] = str(value)
    return MappingProxyType(frozen)


def _freeze_parameters(values: Mapping[str, Any] | None) -> Mapping[str, Any]:
    frozen: dict[str, Any] = {}
    for key, value in (values or {}).items():
        frozen[_validate_text(str(key), "provenance parameter key")] = value
    return MappingProxyType(frozen)


def _freeze_ranges(values: Iterable[tuple[int, int]] | None) -> tuple[tuple[int, int], ...]:
    ranges: list[tuple[int, int]] = []
    for start, end in values or ():
        ranges.append(
            (
                _validate_coordinate(int(start), "affected range start"),
                _validate_coordinate(int(end), "affected range end"),
            )
        )
    return tuple(ranges)


def _utc_timestamp() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")


def _coerce_molecule_type(molecule_type: MoleculeType | str) -> MoleculeType:
    try:
        return molecule_type if isinstance(molecule_type, MoleculeType) else MoleculeType(molecule_type)
    except ValueError as error:
        msg = f"unsupported molecule type: {molecule_type!r}"
        raise ValueError(msg) from error


def _coerce_topology(topology: MoleculeTopology | str) -> MoleculeTopology:
    try:
        return topology if isinstance(topology, MoleculeTopology) else MoleculeTopology(topology)
    except ValueError as error:
        msg = f"unsupported molecule topology: {topology!r}"
        raise ValueError(msg) from error


def _normalize_sequence(sequence: str, molecule_type: MoleculeType) -> str:
    if not isinstance(sequence, str):
        msg = "sequence must be a string"
        raise TypeError(msg)

    normalized = sequence.upper()
    alphabet = {
        MoleculeType.DNA: DNA_ALPHABET,
        MoleculeType.RNA: RNA_ALPHABET,
        MoleculeType.PROTEIN: PROTEIN_ALPHABET,
    }[molecule_type]
    invalid = sorted(set(normalized) - alphabet)
    if invalid:
        invalid_text = "".join(invalid)
        msg = f"invalid {molecule_type.value} sequence characters: {invalid_text}"
        raise ValueError(msg)
    return normalized


def _complement_sequence(sequence: str, molecule_type: MoleculeType) -> str:
    if molecule_type is MoleculeType.DNA:
        return sequence.translate(DNA_COMPLEMENT)
    if molecule_type is MoleculeType.RNA:
        return sequence.translate(RNA_COMPLEMENT)
    msg = "protein records do not support reverse complement"
    raise ValueError(msg)


def _validate_coordinate(value: int, field_name: str) -> int:
    if not isinstance(value, int):
        msg = f"{field_name} must be an integer"
        raise TypeError(msg)
    if value < 0:
        msg = f"{field_name} must be non-negative"
        raise ValueError(msg)
    return value


@dataclass(frozen=True, slots=True)
class FeatureSegment:
    """A zero-based, half-open interval with strand."""

    start: int
    end: int
    strand: int = 0

    def __post_init__(self) -> None:
        start = _validate_coordinate(self.start, "segment start")
        end = _validate_coordinate(self.end, "segment end")
        if end <= start:
            msg = "segment end must be greater than start"
            raise ValueError(msg)
        object.__setattr__(self, "start", start)
        object.__setattr__(self, "end", end)
        object.__setattr__(self, "strand", _validate_strand(self.strand))

    @property
    def length(self) -> int:
        """Return the segment length."""

        return self.end - self.start

    def shift(self, offset: int) -> FeatureSegment:
        """Return this segment shifted by an integer offset."""

        return FeatureSegment(self.start + offset, self.end + offset, self.strand)

    def reverse_complement(self, sequence_length: int) -> FeatureSegment:
        """Map this segment onto the reverse-complemented coordinate system."""

        return FeatureSegment(
            sequence_length - self.end,
            sequence_length - self.start,
            -self.strand if self.strand else 0,
        )


@dataclass(frozen=True, slots=True, init=False)
class Feature:
    """A named annotation that may contain one or more feature segments."""

    type: str
    segments: tuple[FeatureSegment, ...]
    name: str | None = None
    qualifiers: Mapping[str, QualifierValue] = field(default_factory=dict)

    def __init__(
        self,
        *,
        type: str,
        segments: Iterable[FeatureSegment] | None = None,
        start: int | None = None,
        end: int | None = None,
        strand: int = 0,
        name: str | None = None,
        label: str | None = None,
        qualifiers: Mapping[str, Any] | None = None,
    ) -> None:
        if name is not None and label is not None and name != label:
            msg = "feature name and label must match when both are provided"
            raise ValueError(msg)
        feature_name = name if name is not None else label

        if segments is None:
            if start is None or end is None:
                msg = "feature requires either segments or start/end"
                raise ValueError(msg)
            segment_tuple = (FeatureSegment(start, end, strand),)
        else:
            if start is not None or end is not None:
                msg = "feature cannot mix segments with start/end"
                raise ValueError(msg)
            segment_tuple = tuple(segments)
            if not all(isinstance(segment, FeatureSegment) for segment in segment_tuple):
                msg = "feature segments must be FeatureSegment instances"
                raise TypeError(msg)

        if not segment_tuple:
            msg = "feature must contain at least one segment"
            raise ValueError(msg)

        object.__setattr__(self, "type", _validate_text(type, "feature type"))
        object.__setattr__(self, "segments", segment_tuple)
        object.__setattr__(self, "name", _validate_optional_text(feature_name, "feature name"))
        object.__setattr__(self, "qualifiers", _freeze_qualifier_mapping(qualifiers))

    @property
    def label(self) -> str | None:
        """Backward-compatible alias for the feature name."""

        return self.name

    @property
    def start(self) -> int:
        """Return the start of the first segment in feature order."""

        return self.segments[0].start

    @property
    def end(self) -> int:
        """Return the end of the last segment in feature order."""

        return self.segments[-1].end

    @property
    def strand(self) -> int:
        """Return the shared strand, or 0 when segments disagree."""

        strands = {segment.strand for segment in self.segments}
        return strands.pop() if len(strands) == 1 else 0

    def with_segments(self, segments: Iterable[FeatureSegment]) -> Feature:
        """Return a copy of this feature with replacement segments."""

        return Feature(
            type=self.type,
            name=self.name,
            segments=tuple(segments),
            qualifiers=self.qualifiers,
        )

    def validate_for_length(self, sequence_length: int) -> None:
        """Validate this feature for a linear sequence of the given length."""

        self.validate_for_record(sequence_length, MoleculeTopology.LINEAR)

    def validate_for_record(
        self,
        sequence_length: int,
        topology: MoleculeTopology | str,
    ) -> None:
        """Validate segment bounds and compound feature order for a record."""

        topology = _coerce_topology(topology)
        if sequence_length < 0:
            msg = "sequence length must be non-negative"
            raise ValueError(msg)

        for segment in self.segments:
            if segment.end > sequence_length:
                msg = "feature segment end exceeds sequence length"
                raise ValueError(msg)

        sorted_segments = sorted(self.segments, key=lambda segment: (segment.start, segment.end))
        for previous, current in zip(sorted_segments, sorted_segments[1:], strict=False):
            if previous.end > current.start:
                msg = "feature segments must not overlap"
                raise ValueError(msg)

        descents = sum(
            1
            for previous, current in zip(self.segments, self.segments[1:], strict=False)
            if current.start < previous.start
        )
        reverse_strand_transcript_order = (
            all(segment.strand == -1 for segment in self.segments)
            and all(
                current.start < previous.start
                for previous, current in zip(self.segments, self.segments[1:], strict=False)
            )
        )
        if topology is MoleculeTopology.LINEAR and descents and not reverse_strand_transcript_order:
            msg = "linear compound features must be ordered by ascending coordinates"
            raise ValueError(msg)
        if topology is MoleculeTopology.CIRCULAR and descents > 1:
            msg = "circular compound features may wrap around the origin at most once"
            raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class Primer:
    """A primer sequence and optional binding interval."""

    name: str
    sequence: str
    start: int | None = None
    end: int | None = None
    strand: int = 1
    target_id: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", _validate_text(self.name, "primer name"))
        object.__setattr__(self, "sequence", _normalize_sequence(self.sequence, MoleculeType.DNA))
        object.__setattr__(self, "strand", _validate_strand(self.strand))
        object.__setattr__(self, "target_id", _validate_optional_text(self.target_id, "primer target id"))
        if (self.start is None) != (self.end is None):
            msg = "primer start and end must be provided together"
            raise ValueError(msg)
        if self.start is not None and self.end is not None:
            FeatureSegment(self.start, self.end, self.strand)

    def validate_for_length(self, sequence_length: int) -> None:
        """Validate binding coordinates against a target sequence length."""

        if self.end is not None and self.end > sequence_length:
            msg = "primer end exceeds sequence length"
            raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class EnzymeSite:
    """A restriction enzyme recognition site on a sequence record."""

    enzyme_name: str
    recognition_sequence: str
    start: int
    end: int
    strand: int = 1
    cut_index: int | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "enzyme_name", _validate_text(self.enzyme_name, "enzyme name"))
        object.__setattr__(
            self,
            "recognition_sequence",
            _normalize_sequence(self.recognition_sequence, MoleculeType.DNA),
        )
        FeatureSegment(self.start, self.end, self.strand)
        object.__setattr__(self, "strand", _validate_strand(self.strand))
        if self.cut_index is not None:
            cut_index = _validate_coordinate(self.cut_index, "enzyme cut index")
            if cut_index > len(self.recognition_sequence):
                msg = "enzyme cut index exceeds recognition sequence length"
                raise ValueError(msg)

    def validate_for_length(self, sequence_length: int) -> None:
        """Validate site coordinates against a target sequence length."""

        if self.end > sequence_length:
            msg = "enzyme site end exceeds sequence length"
            raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class ProvenanceEvent:
    """Deterministic provenance metadata for a sequence edit."""

    operation: str
    parameters: Mapping[str, Any] = field(default_factory=dict)
    input_record_id: str | None = None
    timestamp: str | None = None
    input_record_ids: tuple[str, ...] = ()
    output_record_id: str | None = None
    input_version_id: str | None = None
    input_version_ids: tuple[str, ...] = ()
    output_version_id: str | None = None
    affected_ranges: tuple[tuple[int, int], ...] = ()
    description: str = ""

    def __post_init__(self) -> None:
        input_ids = tuple(
            _validate_text(str(record_id), "input record id")
            for record_id in self.input_record_ids
        )
        if self.input_record_id is not None and self.input_record_id not in input_ids:
            input_ids = (self.input_record_id, *input_ids)
        input_record_id = self.input_record_id or (input_ids[0] if input_ids else None)

        input_version_ids = tuple(
            _validate_text(str(version_id), "input version id")
            for version_id in self.input_version_ids
        )
        if self.input_version_id is not None and self.input_version_id not in input_version_ids:
            input_version_ids = (self.input_version_id, *input_version_ids)
        input_version_id = self.input_version_id or (
            input_version_ids[0] if input_version_ids else None
        )

        object.__setattr__(self, "operation", _validate_text(self.operation, "operation"))
        object.__setattr__(self, "parameters", _freeze_parameters(self.parameters))
        object.__setattr__(
            self,
            "input_record_id",
            _validate_optional_text(input_record_id, "input record id"),
        )
        object.__setattr__(self, "timestamp", self.timestamp or _utc_timestamp())
        object.__setattr__(self, "input_record_ids", input_ids)
        object.__setattr__(
            self,
            "output_record_id",
            _validate_optional_text(self.output_record_id, "output record id"),
        )
        object.__setattr__(
            self,
            "input_version_id",
            _validate_optional_text(input_version_id, "input version id"),
        )
        object.__setattr__(self, "input_version_ids", input_version_ids)
        object.__setattr__(
            self,
            "output_version_id",
            _validate_optional_text(self.output_version_id, "output version id"),
        )
        object.__setattr__(self, "affected_ranges", _freeze_ranges(self.affected_ranges))
        object.__setattr__(self, "description", str(self.description))

    def with_output(
        self,
        output_record_id: str,
        *,
        output_version_id: str | None = None,
    ) -> ProvenanceEvent:
        """Return a copy with the output record id populated."""

        return ProvenanceEvent(
            operation=self.operation,
            parameters=self.parameters,
            input_record_id=self.input_record_id,
            timestamp=self.timestamp,
            input_record_ids=self.input_record_ids,
            output_record_id=output_record_id,
            input_version_id=self.input_version_id,
            input_version_ids=self.input_version_ids,
            output_version_id=output_version_id or self.output_version_id,
            affected_ranges=self.affected_ranges,
            description=self.description,
        )


@dataclass(frozen=True, slots=True)
class SequenceRecord:
    """Immutable sequence container with explicit molecule type and topology."""

    id: str
    sequence: str
    molecule_type: MoleculeType | str = MoleculeType.DNA
    topology: MoleculeTopology | str = MoleculeTopology.LINEAR
    name: str | None = None
    description: str | None = None
    version_id: str | None = None
    features: tuple[Feature, ...] = ()
    primers: tuple[Primer, ...] = ()
    enzyme_sites: tuple[EnzymeSite, ...] = ()
    history: tuple[ProvenanceEvent, ...] = ()

    def __post_init__(self) -> None:
        molecule_type = _coerce_molecule_type(self.molecule_type)
        topology = _coerce_topology(self.topology)
        sequence = _normalize_sequence(self.sequence, molecule_type)

        object.__setattr__(self, "id", _validate_text(self.id, "record id"))
        object.__setattr__(self, "sequence", sequence)
        object.__setattr__(self, "molecule_type", molecule_type)
        object.__setattr__(self, "topology", topology)
        object.__setattr__(self, "name", _validate_optional_text(self.name, "record name"))
        object.__setattr__(
            self,
            "version_id",
            _validate_optional_text(self.version_id, "record version id")
            or _version_id_from_history(self.id, self.history),
        )
        object.__setattr__(self, "features", tuple(self.features))
        object.__setattr__(self, "primers", tuple(self.primers))
        object.__setattr__(self, "enzyme_sites", tuple(self.enzyme_sites))
        object.__setattr__(self, "history", tuple(self.history))

        for feature in self.features:
            feature.validate_for_record(len(sequence), topology)
        for primer in self.primers:
            primer.validate_for_length(len(sequence))
        for enzyme_site in self.enzyme_sites:
            enzyme_site.validate_for_length(len(sequence))

    @property
    def length(self) -> int:
        """Return the sequence length."""

        return len(self.sequence)

    @property
    def is_circular(self) -> bool:
        """Return whether the record is circular."""

        return self.topology is MoleculeTopology.CIRCULAR

    def reverse_complement(self) -> SequenceRecord:
        """Return the reverse complement of a DNA or RNA record."""

        reversed_sequence = _complement_sequence(self.sequence, self.molecule_type)[::-1]
        reversed_features = tuple(
            feature.with_segments(
                segment.reverse_complement(self.length) for segment in reversed(feature.segments)
            )
            for feature in self.features
        )
        reversed_primers, primer_bindings_remapped = _reverse_complement_primers(
            self.primers,
            self.length,
        )
        annotation_summary = _annotation_summary(
            features_before=len(self.features),
            features_after=len(reversed_features),
            primer_bindings_remapped=primer_bindings_remapped,
            enzyme_sites_dropped=len(self.enzyme_sites),
            primer_policy="binding_coordinates_reverse_complemented",
        )
        return self._copy(
            sequence=reversed_sequence,
            features=reversed_features,
            primers=reversed_primers,
            enzyme_sites=(),
            event=self._event(
                "reverse_complement",
                {"annotation_semantics": annotation_summary},
                affected_ranges=((0, self.length),),
                description=_edit_description(
                    f"Reverse complemented {self.id}",
                    annotation_summary,
                ),
            ),
        )

    def circularize(self) -> SequenceRecord:
        """Return a circular copy of this record."""

        if self.length == 0:
            msg = "cannot circularize an empty sequence"
            raise ValueError(msg)
        primers, invalidated_bindings = _invalidate_primer_bindings(self.primers)
        annotation_summary = _annotation_summary(
            features_before=len(self.features),
            features_after=len(self.features),
            primer_bindings_invalidated=invalidated_bindings,
            enzyme_sites_dropped=len(self.enzyme_sites),
        )
        return self._copy(
            topology=MoleculeTopology.CIRCULAR,
            primers=primers,
            enzyme_sites=(),
            event=self._event(
                "circularize",
                {"annotation_semantics": annotation_summary},
                affected_ranges=((0, self.length),),
                description=_edit_description(f"Circularized {self.id}", annotation_summary),
            ),
        )

    def linearize(self, cut: int = 0) -> SequenceRecord:
        """Return a linear record, cutting a circular sequence at ``cut``."""

        cut = _validate_coordinate(cut, "linearization cut")
        if cut > self.length:
            msg = "linearization cut exceeds sequence length"
            raise ValueError(msg)
        if self.topology is MoleculeTopology.LINEAR:
            if cut != 0:
                msg = "linear records can only be linearized at cut 0"
                raise ValueError(msg)
            primers, invalidated_bindings = _invalidate_primer_bindings(self.primers)
            annotation_summary = _annotation_summary(
                features_before=len(self.features),
                features_after=len(self.features),
                primer_bindings_invalidated=invalidated_bindings,
                enzyme_sites_dropped=len(self.enzyme_sites),
            )
            return self._copy(
                topology=MoleculeTopology.LINEAR,
                primers=primers,
                enzyme_sites=(),
                event=self._event(
                    "linearize",
                    {"cut": cut, "annotation_semantics": annotation_summary},
                    affected_ranges=((0, self.length),),
                    description=_edit_description(
                        f"Linearized {self.id} at {cut}",
                        annotation_summary,
                    ),
                ),
            )

        normalized_cut = 0 if cut == self.length else cut
        intervals = (
            [(normalized_cut, self.length), (0, normalized_cut)]
            if normalized_cut
            else [(0, self.length)]
        )
        linear_sequence = _sequence_from_intervals(self.sequence, intervals)
        linear_features = _project_features(self.features, intervals)
        primers, invalidated_bindings = _invalidate_primer_bindings(self.primers)
        annotation_summary = _annotation_summary(
            features_before=len(self.features),
            features_after=len(linear_features),
            primer_bindings_invalidated=invalidated_bindings,
            enzyme_sites_dropped=len(self.enzyme_sites),
        )
        return self._copy(
            sequence=linear_sequence,
            topology=MoleculeTopology.LINEAR,
            features=linear_features,
            primers=primers,
            enzyme_sites=(),
            event=self._event(
                "linearize",
                {"cut": cut, "annotation_semantics": annotation_summary},
                affected_ranges=((0, self.length),),
                description=_edit_description(f"Linearized {self.id} at {cut}", annotation_summary),
            ),
        )

    def slice(self, start: int, end: int) -> SequenceRecord:
        """Return a linear slice, allowing circular slices that cross the origin."""

        intervals = _selection_intervals(start, end, self.length, self.topology)
        sliced_sequence = _sequence_from_intervals(self.sequence, intervals)
        sliced_features = _project_features(self.features, intervals)
        primers, invalidated_bindings = _invalidate_primer_bindings(self.primers)
        annotation_summary = _annotation_summary(
            features_before=len(self.features),
            features_after=len(sliced_features),
            primer_bindings_invalidated=invalidated_bindings,
            enzyme_sites_dropped=len(self.enzyme_sites),
        )
        return self._copy(
            sequence=sliced_sequence,
            topology=MoleculeTopology.LINEAR,
            features=sliced_features,
            primers=primers,
            enzyme_sites=(),
            event=self._event(
                "slice",
                {"start": start, "end": end, "annotation_semantics": annotation_summary},
                affected_ranges=tuple(intervals),
                description=_edit_description(f"Sliced {self.id}:{start}-{end}", annotation_summary),
            ),
        )

    def insert(self, index: int, insert_sequence: str) -> SequenceRecord:
        """Return a record with ``insert_sequence`` inserted before ``index``."""

        index = _validate_coordinate(index, "insert index")
        if index > self.length:
            msg = "insert index exceeds sequence length"
            raise ValueError(msg)
        normalized_insert = _normalize_sequence(insert_sequence, self.molecule_type)
        insert_length = len(normalized_insert)
        inserted_sequence = self.sequence[:index] + normalized_insert + self.sequence[index:]
        inserted_features = tuple(
            feature.with_segments(_insert_segments(feature.segments, index, insert_length))
            for feature in self.features
        )
        primers, invalidated_bindings = _invalidate_primer_bindings(self.primers)
        annotation_summary = _annotation_summary(
            features_before=len(self.features),
            features_after=len(inserted_features),
            primer_bindings_invalidated=invalidated_bindings,
            enzyme_sites_dropped=len(self.enzyme_sites),
        )
        return self._copy(
            sequence=inserted_sequence,
            features=inserted_features,
            primers=primers,
            enzyme_sites=(),
            event=self._event(
                "insert",
                {
                    "index": index,
                    "inserted_length": insert_length,
                    "annotation_semantics": annotation_summary,
                },
                affected_ranges=((index, index + insert_length),),
                description=_edit_description(
                    f"Inserted {insert_length} bases into {self.id} at {index}",
                    annotation_summary,
                ),
            ),
        )

    def delete(self, start: int, end: int) -> SequenceRecord:
        """Return a record with the selected range removed."""

        kept_intervals = _deletion_complement_intervals(start, end, self.length, self.topology)
        deleted_sequence = _sequence_from_intervals(self.sequence, kept_intervals)
        deleted_features = _project_features(self.features, kept_intervals)
        primers, invalidated_bindings = _invalidate_primer_bindings(self.primers)
        annotation_summary = _annotation_summary(
            features_before=len(self.features),
            features_after=len(deleted_features),
            primer_bindings_invalidated=invalidated_bindings,
            enzyme_sites_dropped=len(self.enzyme_sites),
        )
        return self._copy(
            sequence=deleted_sequence,
            features=deleted_features,
            primers=primers,
            enzyme_sites=(),
            event=self._event(
                "delete",
                {
                    "start": start,
                    "end": end,
                    "deleted_length": self.length - len(deleted_sequence),
                    "annotation_semantics": annotation_summary,
                },
                affected_ranges=tuple(_selection_intervals(start, end, self.length, self.topology)),
                description=_edit_description(
                    f"Deleted {self.length - len(deleted_sequence)} bases from {self.id}",
                    annotation_summary,
                ),
            ),
        )

    def replace(self, start: int, end: int, replacement_sequence: str) -> SequenceRecord:
        """Return a record with the selected range replaced by a new sequence."""

        replacement = _normalize_sequence(replacement_sequence, self.molecule_type)
        replacement_length = len(replacement)

        if self.topology is MoleculeTopology.CIRCULAR and start > end:
            _validate_circular_coordinates(start, end, self.length)
            msg = (
                "origin-crossing circular replace is ambiguous because it would change "
                "the stored sequence origin; delete/insert with explicit coordinates or "
                "linearize at an explicit cut first"
            )
            raise ValueError(msg)

        start, end = _validate_linear_like_range(start, end, self.length, self.topology)
        replaced_sequence = self.sequence[:start] + replacement + self.sequence[end:]
        replaced_features = _replace_features_linear(
            self.features,
            start,
            end,
            replacement_length,
        )
        deleted_length = end - start

        primers, invalidated_bindings = _invalidate_primer_bindings(self.primers)
        annotation_summary = _annotation_summary(
            features_before=len(self.features),
            features_after=len(replaced_features),
            primer_bindings_invalidated=invalidated_bindings,
            enzyme_sites_dropped=len(self.enzyme_sites),
        )
        return self._copy(
            sequence=replaced_sequence,
            features=replaced_features,
            primers=primers,
            enzyme_sites=(),
            event=self._event(
                "replace",
                {
                    "start": start,
                    "end": end,
                    "deleted_length": deleted_length,
                    "inserted_length": replacement_length,
                    "annotation_semantics": annotation_summary,
                },
                affected_ranges=tuple(_selection_intervals(start, end, self.length, self.topology)),
                description=_edit_description(
                    f"Replaced {deleted_length} bases with {replacement_length} bases in {self.id}",
                    annotation_summary,
                ),
            ),
        )

    def _event(
        self,
        operation: str,
        parameters: Mapping[str, Any] | None = None,
        *,
        affected_ranges: Iterable[tuple[int, int]] | None = None,
        description: str = "",
    ) -> ProvenanceEvent:
        return ProvenanceEvent(
            operation=operation,
            parameters=parameters or {},
            input_record_id=self.id,
            input_record_ids=(self.id,),
            input_version_id=self.version_id,
            input_version_ids=(self.version_id,),
            output_record_id=self.id,
            affected_ranges=tuple(affected_ranges or ()),
            description=description,
        )

    def _copy(
        self,
        *,
        sequence: str | None = None,
        topology: MoleculeTopology | str | None = None,
        features: Iterable[Feature] | None = None,
        primers: Iterable[Primer] | None = None,
        enzyme_sites: Iterable[EnzymeSite] | None = None,
        event: ProvenanceEvent | None = None,
    ) -> SequenceRecord:
        final_sequence = self.sequence if sequence is None else sequence
        final_topology = self.topology if topology is None else _coerce_topology(topology)
        final_features = self.features if features is None else tuple(features)
        final_primers = self.primers if primers is None else tuple(primers)
        final_enzyme_sites = self.enzyme_sites if enzyme_sites is None else tuple(enzyme_sites)
        final_version_id = self.version_id
        if event is not None:
            final_version_id = event.output_version_id or _derive_version_id(
                record_id=self.id,
                input_version_id=self.version_id,
                event=event,
                sequence=final_sequence,
                topology=final_topology,
            )
            event = event.with_output(self.id, output_version_id=final_version_id)
        history = self.history + ((event,) if event is not None else ())
        return SequenceRecord(
            id=self.id,
            sequence=final_sequence,
            molecule_type=self.molecule_type,
            topology=final_topology,
            name=self.name,
            description=self.description,
            version_id=final_version_id,
            features=final_features,
            primers=final_primers,
            enzyme_sites=final_enzyme_sites,
            history=history,
        )


def _initial_version_id(record_id: str) -> str:
    return f"{record_id}:v0"


def _version_id_from_history(
    record_id: str,
    history: Iterable[ProvenanceEvent],
) -> str:
    event_tuple = tuple(history)
    for event in reversed(event_tuple):
        if event.output_version_id:
            return event.output_version_id
    return _initial_version_id(record_id)


def _derive_version_id(
    *,
    record_id: str,
    input_version_id: str,
    event: ProvenanceEvent,
    sequence: str,
    topology: MoleculeTopology,
) -> str:
    payload = {
        "record_id": record_id,
        "input_version_id": input_version_id,
        "operation": event.operation,
        "parameters": _json_safe(event.parameters),
        "affected_ranges": event.affected_ranges,
        "sequence": sequence,
        "topology": topology.value,
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()[:16]
    return f"{record_id}:v{digest}"


def _json_safe(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, tuple | list):
        return [_json_safe(item) for item in value]
    if isinstance(value, MoleculeTopology | MoleculeType):
        return value.value
    return value


def _invalidate_primer_bindings(primers: Iterable[Primer]) -> tuple[tuple[Primer, ...], int]:
    """Return primers with any target binding coordinates invalidated."""

    invalidated = 0
    updated: list[Primer] = []
    for primer in primers:
        if primer.start is None:
            updated.append(primer)
            continue
        invalidated += 1
        updated.append(
            Primer(
                name=primer.name,
                sequence=primer.sequence,
                strand=0,
                target_id=primer.target_id,
            )
        )
    return tuple(updated), invalidated


def _reverse_complement_primers(
    primers: Iterable[Primer],
    sequence_length: int,
) -> tuple[tuple[Primer, ...], int]:
    """Map primer binding coordinates onto a reverse-complemented target."""

    remapped = 0
    updated: list[Primer] = []
    for primer in primers:
        if primer.start is None or primer.end is None:
            updated.append(primer)
            continue
        remapped += 1
        updated.append(
            Primer(
                name=primer.name,
                sequence=primer.sequence,
                start=sequence_length - primer.end,
                end=sequence_length - primer.start,
                strand=-primer.strand if primer.strand else 0,
                target_id=primer.target_id,
            )
        )
    return tuple(updated), remapped


def _annotation_summary(
    *,
    features_before: int,
    features_after: int,
    primer_bindings_invalidated: int = 0,
    primer_bindings_remapped: int = 0,
    enzyme_sites_dropped: int = 0,
    primer_policy: str = "binding_coordinates_invalidated",
) -> dict[str, Any]:
    """Return structured provenance for annotation handling during an edit."""

    return {
        "features_before": features_before,
        "features_after": features_after,
        "features_dropped": max(0, features_before - features_after),
        "feature_policy": "remapped_truncated_or_dropped_by_coordinate_projection",
        "primer_policy": primer_policy,
        "primer_bindings_invalidated": primer_bindings_invalidated,
        "primer_bindings_remapped": primer_bindings_remapped,
        "enzyme_site_policy": "dropped_recompute_from_selected_enzyme_set",
        "enzyme_sites_dropped": enzyme_sites_dropped,
    }


def _edit_description(description: str, annotation_summary: Mapping[str, Any]) -> str:
    """Append a concise annotation-handling statement to an edit description."""

    parts = [description, "features remapped/truncated/dropped by coordinate rules"]
    invalidated = int(annotation_summary["primer_bindings_invalidated"])
    remapped = int(annotation_summary["primer_bindings_remapped"])
    dropped_sites = int(annotation_summary["enzyme_sites_dropped"])
    if remapped:
        parts.append(f"{remapped} primer binding(s) reverse-complemented")
    if invalidated:
        parts.append(f"{invalidated} primer binding(s) invalidated for recomputation")
    if dropped_sites:
        parts.append(f"{dropped_sites} enzyme site annotation(s) dropped for recomputation")
    if not remapped and not invalidated:
        parts.append("primer binding annotations unchanged or absent")
    if not dropped_sites:
        parts.append("enzyme site annotations absent")
    return "; ".join(parts)


def _validate_linear_like_range(
    start: int,
    end: int,
    sequence_length: int,
    topology: MoleculeTopology,
) -> tuple[int, int]:
    start = _validate_coordinate(start, "range start")
    end = _validate_coordinate(end, "range end")
    if start > sequence_length or end > sequence_length:
        msg = "range exceeds sequence length"
        raise ValueError(msg)
    if start > end:
        if topology is MoleculeTopology.CIRCULAR:
            msg = "circular ranges that cross the origin must be handled explicitly"
        else:
            msg = "linear range start must be less than or equal to end"
        raise ValueError(msg)
    return start, end


def _validate_circular_coordinates(start: int, end: int, sequence_length: int) -> tuple[int, int]:
    start = _validate_coordinate(start, "range start")
    end = _validate_coordinate(end, "range end")
    if start > sequence_length or end > sequence_length:
        msg = "range exceeds sequence length"
        raise ValueError(msg)
    if sequence_length == 0 and (start != 0 or end != 0):
        msg = "empty circular sequences only support coordinate 0"
        raise ValueError(msg)
    return start, end


def _selection_intervals(
    start: int,
    end: int,
    sequence_length: int,
    topology: MoleculeTopology,
) -> list[tuple[int, int]]:
    if topology is MoleculeTopology.LINEAR:
        start, end = _validate_linear_like_range(start, end, sequence_length, topology)
        return [] if start == end else [(start, end)]

    start, end = _validate_circular_coordinates(start, end, sequence_length)
    if sequence_length == 0:
        return []
    normalized_start = 0 if start == sequence_length else start
    normalized_end = end
    if normalized_start == normalized_end:
        return []
    if normalized_start < normalized_end:
        return [(normalized_start, normalized_end)]
    return [(normalized_start, sequence_length), (0, normalized_end)]


def _deletion_complement_intervals(
    start: int,
    end: int,
    sequence_length: int,
    topology: MoleculeTopology,
) -> list[tuple[int, int]]:
    if topology is MoleculeTopology.LINEAR:
        start, end = _validate_linear_like_range(start, end, sequence_length, topology)
        return _non_empty_intervals((0, start), (end, sequence_length))

    start, end = _validate_circular_coordinates(start, end, sequence_length)
    if sequence_length == 0:
        return []
    normalized_start = 0 if start == sequence_length else start
    normalized_end = end
    if normalized_start == normalized_end:
        return [(0, sequence_length)]
    if normalized_start < normalized_end:
        return _non_empty_intervals((0, normalized_start), (normalized_end, sequence_length))
    return _non_empty_intervals((normalized_end, normalized_start))


def _non_empty_intervals(*intervals: tuple[int, int]) -> list[tuple[int, int]]:
    return [(start, end) for start, end in intervals if start < end]


def _sequence_from_intervals(sequence: str, intervals: Iterable[tuple[int, int]]) -> str:
    return "".join(sequence[start:end] for start, end in intervals)


def _project_features(
    features: Iterable[Feature],
    intervals: Iterable[tuple[int, int]],
    *,
    initial_offset: int = 0,
) -> tuple[Feature, ...]:
    interval_tuple = tuple(intervals)
    projected_features: list[Feature] = []
    for feature in features:
        projected_segments: list[FeatureSegment] = []
        offset = initial_offset
        for interval_start, interval_end in interval_tuple:
            for segment in feature.segments:
                overlap_start = max(segment.start, interval_start)
                overlap_end = min(segment.end, interval_end)
                if overlap_start < overlap_end:
                    projected_segments.append(
                        FeatureSegment(
                            offset + overlap_start - interval_start,
                            offset + overlap_end - interval_start,
                            segment.strand,
                        )
                    )
            offset += interval_end - interval_start
        merged_segments = _merge_sorted_segments(projected_segments)
        if merged_segments:
            projected_features.append(feature.with_segments(merged_segments))
    return tuple(projected_features)


def _merge_sorted_segments(segments: Iterable[FeatureSegment]) -> tuple[FeatureSegment, ...]:
    sorted_segments = sorted(segments, key=lambda segment: (segment.start, segment.end, segment.strand))
    if not sorted_segments:
        return ()

    merged: list[FeatureSegment] = [sorted_segments[0]]
    for segment in sorted_segments[1:]:
        previous = merged[-1]
        if previous.strand == segment.strand and previous.end >= segment.start:
            merged[-1] = FeatureSegment(previous.start, max(previous.end, segment.end), previous.strand)
        else:
            merged.append(segment)
    return tuple(merged)


def _insert_segments(
    segments: Iterable[FeatureSegment],
    index: int,
    insert_length: int,
) -> tuple[FeatureSegment, ...]:
    shifted: list[FeatureSegment] = []
    for segment in segments:
        if segment.end <= index:
            shifted.append(segment)
        elif segment.start >= index:
            shifted.append(FeatureSegment(segment.start + insert_length, segment.end + insert_length, segment.strand))
        else:
            shifted.append(FeatureSegment(segment.start, segment.end + insert_length, segment.strand))
    return tuple(shifted)


def _replace_features_linear(
    features: Iterable[Feature],
    start: int,
    end: int,
    replacement_length: int,
) -> tuple[Feature, ...]:
    deleted_length = end - start
    delta = replacement_length - deleted_length
    replaced_features: list[Feature] = []

    for feature in features:
        replaced_segments: list[FeatureSegment] = []
        for segment in feature.segments:
            if segment.end <= start:
                replaced_segments.append(segment)
            elif segment.start >= end:
                replaced_segments.append(segment.shift(delta))
            else:
                new_start = segment.start if segment.start < start else start
                new_end = segment.end + delta if segment.end > end else start + replacement_length
                if new_start < new_end:
                    replaced_segments.append(FeatureSegment(new_start, new_end, segment.strand))
        merged_segments = _merge_sorted_segments(replaced_segments)
        if merged_segments:
            replaced_features.append(feature.with_segments(merged_segments))
    return tuple(replaced_features)
