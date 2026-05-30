"""Approximate agarose gel lane and band simulation."""

from __future__ import annotations

import json
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from enum import StrEnum
from importlib import resources
from pathlib import Path
from typing import Any


class BandAmountMode(StrEnum):
    """How unspecified fragment amounts are interpreted."""

    EQUAL_MASS = "equal_mass"
    EQUAL_MOLES = "equal_moles"


@dataclass(frozen=True, slots=True)
class GelInputFragment:
    """A DNA fragment loaded onto a gel."""

    size_bp: int
    name: str | None = None
    mass_ng: float | None = None
    moles_fmol: float | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.size_bp, int) or self.size_bp <= 0:
            msg = "fragment size_bp must be a positive integer"
            raise ValueError(msg)
        if self.mass_ng is not None and self.mass_ng < 0:
            msg = "fragment mass_ng must be non-negative"
            raise ValueError(msg)
        if self.moles_fmol is not None and self.moles_fmol < 0:
            msg = "fragment moles_fmol must be non-negative"
            raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class GelLaneInput:
    """Input fragments for one gel lane."""

    name: str
    fragments: tuple[GelInputFragment, ...]

    def __post_init__(self) -> None:
        if not self.name or self.name != self.name.strip():
            msg = "lane name must be non-empty without leading or trailing whitespace"
            raise ValueError(msg)
        object.__setattr__(self, "fragments", tuple(_coerce_fragment(fragment) for fragment in self.fragments))


@dataclass(frozen=True, slots=True)
class LadderBandDefinition:
    """A built-in or user-provided ladder band."""

    size_bp: int
    mass_ng: float | None = None
    label: str | None = None

    def __post_init__(self) -> None:
        if self.size_bp <= 0:
            msg = "ladder band size_bp must be positive"
            raise ValueError(msg)
        if self.mass_ng is not None and self.mass_ng < 0:
            msg = "ladder band mass_ng must be non-negative"
            raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class LadderDefinition:
    """Open JSON ladder definition."""

    id: str
    name: str
    bands: tuple[LadderBandDefinition, ...]
    description: str = ""
    source: str = "PlasmidLab synthetic open definition"

    def __post_init__(self) -> None:
        if not self.id or self.id != self.id.strip():
            msg = "ladder id must be non-empty without leading or trailing whitespace"
            raise ValueError(msg)
        if not self.name or self.name != self.name.strip():
            msg = "ladder name must be non-empty without leading or trailing whitespace"
            raise ValueError(msg)
        if not self.bands:
            msg = "ladder must contain at least one band"
            raise ValueError(msg)
        object.__setattr__(self, "bands", tuple(sorted(self.bands, key=lambda band: band.size_bp, reverse=True)))


@dataclass(frozen=True, slots=True)
class GelBand:
    """A simulated visible gel band."""

    lane_name: str
    size_bp: int
    label: str
    migration: float
    mass_ng: float
    relative_intensity: float
    source_name: str | None = None


@dataclass(frozen=True, slots=True)
class GelLane:
    """Simulated lane output."""

    name: str
    bands: tuple[GelBand, ...]
    is_ladder: bool = False


@dataclass(frozen=True, slots=True)
class AgaroseGelModel:
    """Approximate gel simulation result."""

    lanes: tuple[GelLane, ...]
    agarose_percentage: float
    run_time_minutes: float
    voltage: float
    amount_mode: BandAmountMode
    ladder: LadderDefinition | None = None


def builtin_ladders() -> tuple[LadderDefinition, ...]:
    """Load bundled open JSON ladder definitions."""

    ladder_dir = resources.files("plasmidlab.resources.ladders")
    ladders = []
    for path in sorted(ladder_dir.iterdir(), key=lambda item: item.name):
        if path.suffix == ".json":
            ladders.append(ladder_from_json(path.read_text(encoding="utf-8")))
    return tuple(ladders)


def load_ladder(ladder: str | Path | Mapping[str, Any] | LadderDefinition) -> LadderDefinition:
    """Resolve a built-in, JSON path, mapping, or existing ladder definition."""

    if isinstance(ladder, LadderDefinition):
        return ladder
    if isinstance(ladder, Mapping):
        return ladder_from_mapping(ladder)
    if isinstance(ladder, Path) or (isinstance(ladder, str) and Path(ladder).suffix.lower() == ".json"):
        return ladder_from_json(Path(ladder).read_text(encoding="utf-8"))
    key = str(ladder).strip().lower().replace("_", " ").replace("-", " ")
    for definition in builtin_ladders():
        aliases = {
            definition.id.lower().replace("_", " "),
            definition.name.lower().replace("_", " "),
        }
        if key in aliases:
            return definition
    msg = f"unknown ladder definition: {ladder}"
    raise ValueError(msg)


def ladder_from_json(text: str) -> LadderDefinition:
    """Parse a ladder definition from JSON text."""

    data = json.loads(text)
    if not isinstance(data, Mapping):
        msg = "ladder JSON must contain an object"
        raise ValueError(msg)
    return ladder_from_mapping(data)


def ladder_from_mapping(data: Mapping[str, Any]) -> LadderDefinition:
    """Create a ladder definition from a mapping."""

    bands = tuple(
        LadderBandDefinition(
            size_bp=int(item["size_bp"]),
            mass_ng=None if item.get("mass_ng") is None else float(item["mass_ng"]),
            label=item.get("label"),
        )
        for item in data.get("bands", ())
    )
    return LadderDefinition(
        id=str(data["id"]),
        name=str(data["name"]),
        description=str(data.get("description", "")),
        source=str(data.get("source", "PlasmidLab synthetic open definition")),
        bands=bands,
    )


def simulate_gel(
    lanes: Iterable[GelLaneInput | tuple[str, Iterable[object]] | Iterable[object]],
    *,
    ladder: str | Path | Mapping[str, Any] | LadderDefinition | None = None,
    agarose_percentage: float = 1.0,
    run_time_minutes: float = 45.0,
    voltage: float = 100.0,
    amount_mode: BandAmountMode | str = BandAmountMode.EQUAL_MASS,
    default_band_mass_ng: float = 50.0,
    default_band_moles_fmol: float = 0.05,
) -> AgaroseGelModel:
    """Simulate an approximate agarose gel lane model.

    This is a display model, not an empirical instrument-grade predictor.
    """

    if agarose_percentage <= 0:
        msg = "agarose_percentage must be positive"
        raise ValueError(msg)
    if run_time_minutes <= 0:
        msg = "run_time_minutes must be positive"
        raise ValueError(msg)
    if voltage <= 0:
        msg = "voltage must be positive"
        raise ValueError(msg)
    amount_mode = _coerce_amount_mode(amount_mode)
    lane_inputs = tuple(_coerce_lane_input(index, lane) for index, lane in enumerate(lanes, start=1))
    ladder_definition = load_ladder(ladder) if ladder is not None else None
    if ladder_definition is not None:
        ladder_lane = GelLaneInput(
            name=ladder_definition.name,
            fragments=tuple(
                GelInputFragment(
                    size_bp=band.size_bp,
                    name=band.label,
                    mass_ng=band.mass_ng,
                )
                for band in ladder_definition.bands
            ),
        )
        lane_inputs = (ladder_lane, *lane_inputs)
    if not lane_inputs:
        msg = "gel simulation requires at least one lane or ladder"
        raise ValueError(msg)

    all_fragments = tuple(fragment for lane in lane_inputs for fragment in lane.fragments)
    min_size = min(fragment.size_bp for fragment in all_fragments)
    max_size = max(fragment.size_bp for fragment in all_fragments)
    band_masses = {
        id(fragment): _fragment_mass_ng(
            fragment,
            amount_mode,
            default_band_mass_ng,
            default_band_moles_fmol,
        )
        for fragment in all_fragments
    }
    max_mass = max(band_masses.values(), default=1.0) or 1.0

    output_lanes: list[GelLane] = []
    for lane_index, lane in enumerate(lane_inputs):
        bands = tuple(
            GelBand(
                lane_name=lane.name,
                size_bp=fragment.size_bp,
                label=fragment.name or _format_size(fragment.size_bp),
                migration=_migration(
                    fragment.size_bp,
                    min_size,
                    max_size,
                    agarose_percentage,
                    run_time_minutes,
                    voltage,
                ),
                mass_ng=band_masses[id(fragment)],
                relative_intensity=band_masses[id(fragment)] / max_mass,
                source_name=fragment.name,
            )
            for fragment in sorted(lane.fragments, key=lambda item: item.size_bp, reverse=True)
        )
        output_lanes.append(
            GelLane(
                name=lane.name,
                bands=bands,
                is_ladder=ladder_definition is not None and lane_index == 0,
            )
        )

    return AgaroseGelModel(
        lanes=tuple(output_lanes),
        agarose_percentage=float(agarose_percentage),
        run_time_minutes=float(run_time_minutes),
        voltage=float(voltage),
        amount_mode=amount_mode,
        ladder=ladder_definition,
    )


def _coerce_lane_input(
    index: int,
    lane: GelLaneInput | tuple[str, Iterable[object]] | Iterable[object],
) -> GelLaneInput:
    if isinstance(lane, GelLaneInput):
        return lane
    if isinstance(lane, tuple) and len(lane) == 2 and isinstance(lane[0], str):
        return GelLaneInput(name=lane[0], fragments=tuple(_coerce_fragment(fragment) for fragment in lane[1]))
    return GelLaneInput(name=f"Lane {index}", fragments=tuple(_coerce_fragment(fragment) for fragment in lane))


def _coerce_fragment(fragment: object) -> GelInputFragment:
    if isinstance(fragment, GelInputFragment):
        return fragment
    if isinstance(fragment, int):
        return GelInputFragment(size_bp=fragment)
    size = getattr(fragment, "length", None)
    if isinstance(size, int):
        source_names = getattr(fragment, "source_enzymes", ())
        name = ",".join(source_names) if source_names else None
        return GelInputFragment(size_bp=size, name=name)
    sequence = getattr(fragment, "sequence", None)
    if isinstance(sequence, str):
        return GelInputFragment(size_bp=len(sequence))
    msg = f"unsupported gel fragment input: {fragment!r}"
    raise TypeError(msg)


def _coerce_amount_mode(mode: BandAmountMode | str) -> BandAmountMode:
    try:
        return mode if isinstance(mode, BandAmountMode) else BandAmountMode(mode)
    except ValueError as error:
        msg = f"unsupported gel amount mode: {mode!r}"
        raise ValueError(msg) from error


def _fragment_mass_ng(
    fragment: GelInputFragment,
    amount_mode: BandAmountMode,
    default_band_mass_ng: float,
    default_band_moles_fmol: float,
) -> float:
    if fragment.mass_ng is not None:
        return float(fragment.mass_ng)
    if fragment.moles_fmol is not None:
        return _fmol_to_ng(fragment.moles_fmol, fragment.size_bp)
    if amount_mode is BandAmountMode.EQUAL_MASS:
        return float(default_band_mass_ng)
    return _fmol_to_ng(default_band_moles_fmol, fragment.size_bp)


def _fmol_to_ng(fmol: float, size_bp: int) -> float:
    return float(fmol) * float(size_bp) * 0.00066


def _migration(
    size_bp: int,
    min_size: int,
    max_size: int,
    agarose_percentage: float,
    run_time_minutes: float,
    voltage: float,
) -> float:
    if min_size == max_size:
        base = 0.5
    else:
        import math

        log_min = math.log10(min_size)
        log_max = math.log10(max_size)
        base = (log_max - math.log10(size_bp)) / (log_max - log_min)
    agarose_exponent = 0.9 + max(0.0, agarose_percentage - 1.0) * 0.16
    run_scale = min(1.2, max(0.45, (run_time_minutes / 45.0 * voltage / 100.0) ** 0.5))
    migration = (base**agarose_exponent) * run_scale
    return min(0.98, max(0.02, migration))


def _format_size(size_bp: int) -> str:
    if size_bp >= 1000:
        value = size_bp / 1000
        if value.is_integer():
            return f"{int(value)} kb"
        return f"{value:.1f} kb"
    return f"{size_bp} bp"
