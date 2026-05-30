"""Automatic feature detection from open, user-editable libraries."""

from __future__ import annotations

import json
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from enum import StrEnum
from importlib import resources
from pathlib import Path
from typing import Any

from plasmidlab.core.models import (
    Feature,
    FeatureSegment,
    MoleculeTopology,
    MoleculeType,
    ProvenanceEvent,
    SequenceRecord,
)


IUPAC_DNA: Mapping[str, frozenset[str]] = {
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
DEFAULT_FEATURE_LIBRARY_PACKAGE = "plasmidlab.data.features"


class FeatureStrandBehavior(StrEnum):
    """How a library entry should be searched and annotated."""

    BOTH = "both"
    FORWARD_ONLY = "forward_only"
    REVERSE_ONLY = "reverse_only"
    NON_DIRECTIONAL = "non_directional"


@dataclass(frozen=True, slots=True)
class FeatureLibraryEntry:
    """One open feature-library entry."""

    name: str
    type: str
    sequence: str
    aliases: tuple[str, ...] = ()
    strand_behavior: FeatureStrandBehavior = FeatureStrandBehavior.BOTH
    minimum_identity: float = 1.0
    notes: str = ""
    source_path: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", _validate_text(self.name, "feature name"))
        object.__setattr__(self, "type", _validate_text(self.type, "feature type"))
        object.__setattr__(self, "sequence", _normalize_dna(self.sequence, "feature sequence"))
        object.__setattr__(
            self,
            "aliases",
            tuple(_validate_text(str(alias), "feature alias") for alias in self.aliases),
        )
        object.__setattr__(self, "strand_behavior", _coerce_strand_behavior(self.strand_behavior))
        identity = float(self.minimum_identity)
        if not 0 < identity <= 1:
            msg = "minimum_identity must be greater than 0 and less than or equal to 1"
            raise ValueError(msg)
        object.__setattr__(self, "minimum_identity", identity)
        object.__setattr__(self, "notes", str(self.notes))
        if self.source_path is not None:
            object.__setattr__(self, "source_path", str(self.source_path))

    @property
    def length(self) -> int:
        """Return the library sequence length."""

        return len(self.sequence)

    @classmethod
    def from_mapping(
        cls,
        data: Mapping[str, Any],
        *,
        source_path: str | None = None,
    ) -> FeatureLibraryEntry:
        """Create an entry from a JSON object."""

        return cls(
            name=str(data["name"]),
            type=str(data["type"]),
            sequence=str(data["sequence"]),
            aliases=tuple(str(alias) for alias in data.get("aliases", ())),
            strand_behavior=_coerce_strand_behavior(data.get("strand_behavior", FeatureStrandBehavior.BOTH)),
            minimum_identity=float(data.get("minimum_identity", 1.0)),
            notes=str(data.get("notes", "")),
            source_path=source_path,
        )


@dataclass(frozen=True, slots=True)
class FeatureDetection:
    """A proposed annotation from a feature-library match."""

    entry: FeatureLibraryEntry
    start: int
    end: int
    segments: tuple[FeatureSegment, ...]
    strand: int
    identity: float
    mismatches: int
    matched_sequence: str
    wraps_origin: bool = False

    @property
    def name(self) -> str:
        """Return the proposed feature name."""

        return self.entry.name

    @property
    def type(self) -> str:
        """Return the proposed feature type."""

        return self.entry.type

    @property
    def length(self) -> int:
        """Return the matched feature length."""

        return sum(segment.length for segment in self.segments)

    @property
    def exact(self) -> bool:
        """Return whether the match is exact."""

        return self.mismatches == 0

    def to_feature(self) -> Feature:
        """Convert this proposal into a core Feature annotation."""

        qualifiers: dict[str, str] = {
            "label": self.entry.name,
            "detected_identity": f"{self.identity:.3f}",
            "detection": "PlasmidLab open feature library",
        }
        if self.entry.aliases:
            qualifiers["aliases"] = "; ".join(self.entry.aliases)
        if self.entry.notes:
            qualifiers["note"] = self.entry.notes
        return Feature(
            type=self.entry.type,
            name=self.entry.name,
            segments=self.segments,
            qualifiers=qualifiers,
        )


def load_feature_library(
    paths: str | Path | Iterable[str | Path] | None = None,
) -> tuple[FeatureLibraryEntry, ...]:
    """Load feature-library JSON files from files or directories.

    When ``paths`` is omitted, PlasmidLab loads its packaged open default library via
    ``importlib.resources`` so it works from source checkouts, editable installs, wheels,
    and PyInstaller bundles. Pass explicit files or directories to load user-provided
    lab libraries.
    """

    entries: list[FeatureLibraryEntry] = []
    if paths is None:
        for source_name, raw_json in _default_library_json_documents():
            data = json.loads(raw_json)
            entries.extend(_entries_from_json_data(data, source_path=source_name))
        return tuple(entries)

    resolved_paths = _resolve_library_paths(paths)
    for path in resolved_paths:
        data = json.loads(path.read_text(encoding="utf-8"))
        entries.extend(_entries_from_json_data(data, source_path=str(path)))
    return tuple(entries)


def detect_features(
    record: SequenceRecord,
    library: Iterable[FeatureLibraryEntry] | None = None,
    *,
    library_paths: str | Path | Iterable[str | Path] | None = None,
    resolve_duplicates: bool = True,
) -> tuple[FeatureDetection, ...]:
    """Detect exact and near-exact library features in a DNA record."""

    if record.molecule_type is not MoleculeType.DNA:
        msg = "feature detection requires a DNA sequence"
        raise ValueError(msg)
    entries = tuple(library) if library is not None else load_feature_library(library_paths)
    detections = tuple(
        detection
        for entry in entries
        for detection in _detect_entry(record, entry)
    )
    if resolve_duplicates:
        detections = resolve_duplicate_detections(detections)
    return tuple(sorted(detections, key=lambda item: (item.start, item.end, item.name, -item.identity)))


def resolve_duplicate_detections(
    detections: Iterable[FeatureDetection],
) -> tuple[FeatureDetection, ...]:
    """Suppress repeated proposals for the same site, keeping the best match."""

    kept: list[FeatureDetection] = []
    for candidate in sorted(detections, key=lambda item: (-item.identity, -item.length, item.start, item.name)):
        if any(_is_duplicate_detection(candidate, previous) for previous in kept):
            continue
        kept.append(candidate)
    return tuple(sorted(kept, key=lambda item: (item.start, item.end, item.name)))


def apply_feature_annotations(
    record: SequenceRecord,
    detections: Iterable[FeatureDetection],
) -> SequenceRecord:
    """Return a new record with selected feature detections applied."""

    selected = tuple(detections)
    if not selected:
        return record
    existing_keys = {_feature_key(feature) for feature in record.features}
    new_features: list[Feature] = []
    for detection in selected:
        feature = detection.to_feature()
        key = _feature_key(feature)
        if key in existing_keys:
            continue
        new_features.append(feature)
        existing_keys.add(key)
    if not new_features:
        return record

    event = ProvenanceEvent(
        operation="detect_features",
        parameters={
            "applied_count": len(new_features),
            "names": tuple(feature.name or feature.type for feature in new_features),
        },
        input_record_id=record.id,
        input_record_ids=(record.id,),
        input_version_id=record.version_id,
        input_version_ids=(record.version_id,),
        output_record_id=record.id,
        affected_ranges=tuple(
            (segment.start, segment.end)
            for feature in new_features
            for segment in feature.segments
        ),
        description=f"Detected and applied {len(new_features)} features on {record.id}",
    )
    return record._copy(features=record.features + tuple(new_features), event=event)


def _detect_entry(record: SequenceRecord, entry: FeatureLibraryEntry) -> tuple[FeatureDetection, ...]:
    if entry.length == 0 or record.length == 0:
        return ()
    if record.topology is MoleculeTopology.LINEAR and entry.length > record.length:
        return ()
    if record.topology is MoleculeTopology.CIRCULAR and entry.length > record.length:
        return ()

    patterns = _patterns_for_entry(entry)
    max_start = record.length if record.topology is MoleculeTopology.CIRCULAR else record.length - entry.length + 1
    detections: list[FeatureDetection] = []
    for pattern, strand in patterns:
        for start in range(max_start):
            matched_sequence = _template_window(record.sequence, start, entry.length, record.topology)
            if matched_sequence is None:
                continue
            identity, mismatches = _identity(matched_sequence, pattern)
            if identity + 1e-12 < entry.minimum_identity:
                continue
            raw_end = start + entry.length
            wraps_origin = record.topology is MoleculeTopology.CIRCULAR and raw_end > record.length
            end = raw_end % record.length if wraps_origin else raw_end
            detections.append(
                FeatureDetection(
                    entry=entry,
                    start=start,
                    end=end,
                    segments=_segments_for_match(start, entry.length, strand, record.length, wraps_origin),
                    strand=strand,
                    identity=identity,
                    mismatches=mismatches,
                    matched_sequence=matched_sequence,
                    wraps_origin=wraps_origin,
                )
            )
    return tuple(detections)


def _patterns_for_entry(entry: FeatureLibraryEntry) -> tuple[tuple[str, int], ...]:
    if entry.strand_behavior is FeatureStrandBehavior.FORWARD_ONLY:
        return ((entry.sequence, 1),)
    if entry.strand_behavior is FeatureStrandBehavior.REVERSE_ONLY:
        return ((_reverse_complement(entry.sequence), -1),)
    if entry.strand_behavior is FeatureStrandBehavior.NON_DIRECTIONAL:
        reverse = _reverse_complement(entry.sequence)
        if reverse == entry.sequence:
            return ((entry.sequence, 0),)
        return ((entry.sequence, 0), (reverse, 0))
    reverse = _reverse_complement(entry.sequence)
    if reverse == entry.sequence:
        return ((entry.sequence, 1),)
    return ((entry.sequence, 1), (reverse, -1))


def _template_window(sequence: str, start: int, length: int, topology: MoleculeTopology) -> str | None:
    end = start + length
    if topology is MoleculeTopology.LINEAR:
        return sequence[start:end] if end <= len(sequence) else None
    if end <= len(sequence):
        return sequence[start:end]
    return sequence[start:] + sequence[: end % len(sequence)]


def _segments_for_match(
    start: int,
    length: int,
    strand: int,
    record_length: int,
    wraps_origin: bool,
) -> tuple[FeatureSegment, ...]:
    if not wraps_origin:
        return (FeatureSegment(start, start + length, strand),)
    return (
        FeatureSegment(start, record_length, strand),
        FeatureSegment(0, (start + length) % record_length, strand),
    )


def _identity(template: str, pattern: str) -> tuple[float, int]:
    matches = sum(
        1
        for template_base, pattern_base in zip(template, pattern, strict=True)
        if IUPAC_DNA[template_base] & IUPAC_DNA[pattern_base]
    )
    mismatches = len(pattern) - matches
    return matches / len(pattern), mismatches


def _is_duplicate_detection(candidate: FeatureDetection, previous: FeatureDetection) -> bool:
    if _canonical_segments(candidate) == _canonical_segments(previous) and candidate.strand == previous.strand:
        return True
    if candidate.name != previous.name:
        return False
    overlap = len(_covered_positions(candidate) & _covered_positions(previous))
    return overlap / max(1, min(candidate.length, previous.length)) >= 0.8


def _covered_positions(detection: FeatureDetection) -> set[int]:
    positions: set[int] = set()
    for segment in detection.segments:
        positions.update(range(segment.start, segment.end))
    return positions


def _canonical_segments(detection: FeatureDetection) -> tuple[tuple[int, int], ...]:
    return tuple((segment.start, segment.end) for segment in detection.segments)


def _feature_key(feature: Feature) -> tuple[str, str | None, tuple[tuple[int, int, int], ...]]:
    return (
        feature.type,
        feature.name,
        tuple((segment.start, segment.end, segment.strand) for segment in feature.segments),
    )


def _resolve_library_paths(paths: str | Path | Iterable[str | Path] | None) -> tuple[Path, ...]:
    if paths is None:
        return ()
    if isinstance(paths, str | Path):
        roots = (Path(paths),)
    else:
        roots = tuple(Path(path) for path in paths)
    resolved: list[Path] = []
    for root in roots:
        if root.is_dir():
            resolved.extend(sorted(root.glob("*.json"), key=lambda path: path.name))
        elif root.is_file():
            resolved.append(root)
    return tuple(dict.fromkeys(path.resolve() for path in resolved))


def _default_library_json_documents() -> tuple[tuple[str, str], ...]:
    root = resources.files(DEFAULT_FEATURE_LIBRARY_PACKAGE)
    documents: list[tuple[str, str]] = []
    for resource in sorted(root.iterdir(), key=lambda item: item.name):
        if resource.is_file() and resource.name.endswith(".json"):
            source_name = f"package:{DEFAULT_FEATURE_LIBRARY_PACKAGE}/{resource.name}"
            documents.append((source_name, resource.read_text(encoding="utf-8")))
    return tuple(documents)


def _entries_from_json_data(data: Any, *, source_path: str) -> tuple[FeatureLibraryEntry, ...]:
    if isinstance(data, Mapping):
        if "features" in data:
            raw_entries = data["features"]
        else:
            raw_entries = (data,)
    elif isinstance(data, list):
        raw_entries = data
    else:
        msg = "feature library JSON must contain an object, a features list, or a list"
        raise ValueError(msg)
    return tuple(
        FeatureLibraryEntry.from_mapping(_require_mapping(item), source_path=source_path)
        for item in raw_entries
    )


def _require_mapping(value: Any) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        msg = "feature library entries must be JSON objects"
        raise ValueError(msg)
    return value


def _coerce_strand_behavior(value: FeatureStrandBehavior | str) -> FeatureStrandBehavior:
    if isinstance(value, FeatureStrandBehavior):
        return value
    normalized = str(value).strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "both_strands": "both",
        "forward": "forward_only",
        "plus": "forward_only",
        "reverse": "reverse_only",
        "minus": "reverse_only",
        "none": "non_directional",
        "unknown": "non_directional",
    }
    normalized = aliases.get(normalized, normalized)
    try:
        return FeatureStrandBehavior(normalized)
    except ValueError as error:
        msg = f"unsupported strand_behavior: {value!r}"
        raise ValueError(msg) from error


def _validate_text(value: str, field_name: str) -> str:
    if not isinstance(value, str):
        msg = f"{field_name} must be a string"
        raise TypeError(msg)
    if not value or value != value.strip():
        msg = f"{field_name} must be non-empty without leading or trailing whitespace"
        raise ValueError(msg)
    return value


def _normalize_dna(sequence: str, field_name: str) -> str:
    normalized = "".join(str(sequence).upper().split()).replace("U", "T")
    if not normalized:
        msg = f"{field_name} must not be empty"
        raise ValueError(msg)
    invalid = sorted(set(normalized) - set(IUPAC_DNA))
    if invalid:
        msg = f"invalid DNA characters in {field_name}: {''.join(invalid)}"
        raise ValueError(msg)
    return normalized


def _reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(DNA_COMPLEMENT)[::-1]
