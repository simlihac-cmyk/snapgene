"""GenBank import and export for PlasmidLab sequence records."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from io import StringIO
from pathlib import Path
from typing import Any, TextIO
import warnings

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, ExactPosition, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord as BioSeqRecord

from plasmidlab.core import (
    Feature,
    FeatureSegment,
    MoleculeTopology,
    MoleculeType,
    SequenceRecord,
)

LOCATION_OPERATOR_QUALIFIER = "_plasmidlab_location_operator"
LOCATION_WARNINGS_QUALIFIER = "_plasmidlab_location_warnings"
FUZZY_LOCATION_QUALIFIER = "_plasmidlab_fuzzy_location"
INTERNAL_QUALIFIER_PREFIX = "_plasmidlab_"
SUPPORTED_COMPOUND_OPERATORS = frozenset({"join", "order"})


class GenBankLossyWarning(UserWarning):
    """Warning emitted when GenBank data is imported or exported with explicit degradation."""


def loads_genbank(text: str) -> tuple[SequenceRecord, ...]:
    """Parse GenBank text into PlasmidLab records."""

    return parse_genbank(StringIO(text))


def dumps_genbank(records: SequenceRecord | Iterable[SequenceRecord]) -> str:
    """Serialize one or more PlasmidLab records to GenBank text."""

    output = StringIO()
    write_genbank(records, output)
    return output.getvalue()


def read_genbank(path: str | Path) -> tuple[SequenceRecord, ...]:
    """Read GenBank records from a path."""

    with Path(path).open(encoding="utf-8") as handle:
        return parse_genbank(handle)


def write_genbank(records: SequenceRecord | Iterable[SequenceRecord], target: str | Path | TextIO) -> None:
    """Write one or more PlasmidLab records to GenBank."""

    bio_records = [to_biopython_record(record) for record in _as_record_tuple(records)]
    if hasattr(target, "write"):
        SeqIO.write(bio_records, target, "genbank")
        return

    with Path(target).open("w", encoding="utf-8") as handle:
        SeqIO.write(bio_records, handle, "genbank")


def parse_genbank(handle: TextIO) -> tuple[SequenceRecord, ...]:
    """Parse GenBank records from a text handle."""

    return tuple(from_biopython_record(bio_record) for bio_record in SeqIO.parse(handle, "genbank"))


def from_biopython_record(bio_record: BioSeqRecord) -> SequenceRecord:
    """Convert a Biopython SeqRecord into a PlasmidLab SequenceRecord."""

    molecule_type = _molecule_type_from_annotations(bio_record.annotations)
    topology = _topology_from_annotations(bio_record.annotations)
    return SequenceRecord(
        id=bio_record.id,
        name=bio_record.name if bio_record.name and bio_record.name != "<unknown name>" else None,
        description=bio_record.description if bio_record.description != "<unknown description>" else None,
        sequence=str(bio_record.seq),
        molecule_type=molecule_type,
        topology=topology,
        features=tuple(_feature_from_biopython(feature) for feature in bio_record.features),
    )


def to_biopython_record(record: SequenceRecord) -> BioSeqRecord:
    """Convert a PlasmidLab SequenceRecord into a Biopython SeqRecord."""

    bio_record = BioSeqRecord(
        Seq(record.sequence),
        id=record.id,
        name=record.name or record.id,
        description=record.description or record.name or record.id,
    )
    bio_record.annotations["molecule_type"] = _molecule_type_to_genbank(record.molecule_type)
    bio_record.annotations["topology"] = record.topology.value
    bio_record.annotations["data_file_division"] = "SYN"
    bio_record.annotations["date"] = "01-JAN-2000"
    bio_record.features = [_feature_to_biopython(feature) for feature in record.features]
    return bio_record


def _feature_from_biopython(bio_feature: SeqFeature) -> Feature:
    if bio_feature.location is None:
        msg = f"GenBank feature {bio_feature.type!r} has no location"
        raise ValueError(msg)

    _reject_remote_locations(bio_feature.location, bio_feature.type)

    parts = _location_parts(bio_feature.location)
    segments: list[FeatureSegment] = []
    lossy_warnings: list[str] = []
    fuzzy_location = False
    for part in parts:
        if _has_fuzzy_boundary(part):
            fuzzy_location = True
            lossy_warnings.append(
                "GenBank fuzzy location for feature "
                f"{bio_feature.type!r} imported with numeric boundaries: {part}"
            )
        segments.append(
            FeatureSegment(
                start=_numeric_position(part.start, "feature start"),
                end=_numeric_position(part.end, "feature end"),
                strand=_strand_from_biopython(part.strand),
            )
        )

    qualifiers = _qualifiers_from_biopython(bio_feature.qualifiers)
    operator = getattr(bio_feature.location, "operator", None)
    if operator is not None:
        operator = str(operator)
        if operator in SUPPORTED_COMPOUND_OPERATORS:
            qualifiers[LOCATION_OPERATOR_QUALIFIER] = operator
        else:
            lossy_warnings.append(
                f"GenBank compound location operator {operator!r} is not supported for "
                f"feature {bio_feature.type!r}; segment coordinates were imported"
            )

    if lossy_warnings:
        for message in lossy_warnings:
            _warn_lossy(message)
        qualifiers[LOCATION_WARNINGS_QUALIFIER] = (
            *_qualifier_values(qualifiers.get(LOCATION_WARNINGS_QUALIFIER)),
            *lossy_warnings,
        )
        if fuzzy_location:
            qualifiers[FUZZY_LOCATION_QUALIFIER] = str(bio_feature.location)

    return Feature(
        type=bio_feature.type,
        name=_name_from_qualifiers(qualifiers),
        segments=tuple(segments),
        qualifiers=qualifiers,
    )


def _feature_to_biopython(feature: Feature) -> SeqFeature:
    if FUZZY_LOCATION_QUALIFIER in feature.qualifiers:
        _warn_lossy(
            f"Feature {feature.name or feature.type!r} was imported from a fuzzy GenBank "
            "location; export uses exact numeric boundaries because fuzzy boundary metadata "
            "is retained only as PlasmidLab warning metadata"
        )

    parts = [
        SimpleLocation(
            segment.start,
            segment.end,
            strand=_strand_to_biopython(segment.strand),
        )
        for segment in feature.segments
    ]
    if len(parts) == 1:
        location = parts[0]
        if LOCATION_OPERATOR_QUALIFIER in feature.qualifiers:
            _warn_lossy(
                f"Feature {feature.name or feature.type!r} has a compound operator marker "
                "but only one segment; the operator marker was ignored on export"
            )
    else:
        location = CompoundLocation(parts, operator=_compound_operator_for_feature(feature))
    return SeqFeature(
        location=location,
        type=feature.type,
        qualifiers=_qualifiers_to_biopython(feature),
    )


def _location_parts(location: Any) -> tuple[Any, ...]:
    parts = getattr(location, "parts", None)
    if parts is None:
        return (location,)
    return tuple(parts)


def _reject_remote_locations(location: Any, feature_type: str) -> None:
    for part in _location_parts(location):
        ref = getattr(part, "ref", None)
        ref_db = getattr(part, "ref_db", None)
        if ref or ref_db:
            msg = (
                "remote GenBank locations are not supported for feature "
                f"{feature_type!r}: {part}"
            )
            raise ValueError(msg)


def _has_fuzzy_boundary(location: Any) -> bool:
    return not isinstance(location.start, ExactPosition) or not isinstance(location.end, ExactPosition)


def _numeric_position(position: Any, field_name: str) -> int:
    try:
        return int(position)
    except (TypeError, ValueError) as error:
        msg = f"unsupported GenBank {field_name}: {position!r}"
        raise ValueError(msg) from error


def _strand_from_biopython(strand: int | None) -> int:
    return strand if strand in (-1, 1) else 0


def _strand_to_biopython(strand: int) -> int | None:
    return strand if strand in (-1, 1) else None


def _qualifiers_from_biopython(qualifiers: Mapping[str, Any]) -> dict[str, str | tuple[str, ...]]:
    normalized: dict[str, str | tuple[str, ...]] = {}
    for key, value in qualifiers.items():
        if isinstance(value, list):
            normalized[str(key)] = tuple(str(item) for item in value)
        else:
            normalized[str(key)] = str(value)
    return normalized


def _qualifiers_to_biopython(feature: Feature) -> dict[str, list[str]]:
    qualifiers: dict[str, list[str]] = {}
    for key, value in feature.qualifiers.items():
        if key.startswith(INTERNAL_QUALIFIER_PREFIX):
            continue
        qualifiers[key] = list(_qualifier_values(value))
    if feature.name and not any(key in qualifiers for key in ("label", "gene", "product", "note")):
        qualifiers["label"] = [feature.name]
    return qualifiers


def _compound_operator_for_feature(feature: Feature) -> str:
    operator = _first_qualifier_value(feature.qualifiers.get(LOCATION_OPERATOR_QUALIFIER))
    if not operator:
        return "join"
    if operator not in SUPPORTED_COMPOUND_OPERATORS:
        _warn_lossy(
            f"Feature {feature.name or feature.type!r} has unsupported compound location "
            f"operator {operator!r}; exporting as join"
        )
        return "join"
    return operator


def _name_from_qualifiers(qualifiers: Mapping[str, Any]) -> str | None:
    for key in ("label", "gene", "product", "note"):
        for value in _qualifier_values(qualifiers.get(key)):
            stripped = value.strip()
            if stripped:
                return stripped
    return None


def _qualifier_values(value: Any) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        return (value,)
    if isinstance(value, tuple | list):
        return tuple(str(item) for item in value)
    return (str(value),)


def _first_qualifier_value(value: Any) -> str | None:
    values = _qualifier_values(value)
    return values[0] if values else None


def _warn_lossy(message: str) -> None:
    warnings.warn(message, GenBankLossyWarning, stacklevel=3)


def _molecule_type_from_annotations(annotations: Mapping[str, Any]) -> MoleculeType:
    molecule_type = str(annotations.get("molecule_type", "DNA")).lower()
    if "rna" in molecule_type:
        return MoleculeType.RNA
    if "protein" in molecule_type or "amino" in molecule_type:
        return MoleculeType.PROTEIN
    return MoleculeType.DNA


def _molecule_type_to_genbank(molecule_type: MoleculeType) -> str:
    if molecule_type is MoleculeType.PROTEIN:
        return "protein"
    return molecule_type.value


def _topology_from_annotations(annotations: Mapping[str, Any]) -> MoleculeTopology:
    topology = str(annotations.get("topology", "linear")).lower()
    if topology == MoleculeTopology.CIRCULAR.value:
        return MoleculeTopology.CIRCULAR
    return MoleculeTopology.LINEAR


def _as_record_tuple(records: SequenceRecord | Iterable[SequenceRecord]) -> tuple[SequenceRecord, ...]:
    if isinstance(records, SequenceRecord):
        return (records,)
    return tuple(records)
