"""GUI-independent plasmid map drawing model and SVG export."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from math import cos, pi, sin
from pathlib import Path
from typing import Literal
from xml.sax.saxutils import escape

from plasmidlab.core import (
    EnzymeSite,
    Feature,
    FeatureSegment,
    MoleculeTopology,
    Primer,
    PrimerBinding,
    RestrictionAnalysis,
    RestrictionSite,
    SequenceRecord,
)


MapKind = Literal["circular", "linear"]


DEFAULT_FEATURE_COLORS: Mapping[str, str] = {
    "promoter": "#2f9e44",
    "CDS": "#1c7ed6",
    "terminator": "#e03131",
    "rep_origin": "#7048e8",
    "origin": "#7048e8",
    "misc_feature": "#868e96",
    "primer": "#f08c00",
    "enzyme": "#212529",
}


@dataclass(frozen=True, slots=True)
class MapStyle:
    """Layout and style settings for plasmid map rendering."""

    width: int = 900
    height: int = 700
    margin: int = 72
    circular_radius: int = 210
    feature_lane_gap: int = 20
    feature_stroke_width: int = 12
    primer_lane_gap: int = 16
    label_font_size: int = 12
    sequence_stroke: str = "#343a40"
    label_color: str = "#212529"
    background_color: str = "#ffffff"
    feature_colors: Mapping[str, str] = field(default_factory=lambda: dict(DEFAULT_FEATURE_COLORS))


@dataclass(frozen=True, slots=True)
class Point:
    """A 2D drawing point."""

    x: float
    y: float


@dataclass(frozen=True, slots=True)
class FeatureArc:
    """A feature interval drawn as an arc or linear segment."""

    feature_index: int
    segment_index: int
    feature_type: str
    label: str
    start: int
    end: int
    strand: int
    lane: int
    color: str
    map_kind: MapKind
    start_angle: float | None = None
    end_angle: float | None = None
    radius: float | None = None
    center: Point | None = None
    start_point: Point | None = None
    end_point: Point | None = None


@dataclass(frozen=True, slots=True)
class FeatureArrow:
    """Arrowhead for a directional feature."""

    feature_index: int
    segment_index: int
    strand: int
    color: str
    points: tuple[Point, Point, Point]


@dataclass(frozen=True, slots=True)
class EnzymeTick:
    """Restriction enzyme tick mark."""

    enzyme_name: str
    position: int
    recognition_sequence: str
    start_point: Point
    end_point: Point
    angle: float | None = None
    site_start: int | None = None
    site_end: int | None = None
    wraps_origin: bool = False


@dataclass(frozen=True, slots=True)
class PrimerArrow:
    """Primer binding arrow."""

    name: str
    start: int
    end: int
    strand: int
    color: str
    points: tuple[Point, Point, Point]
    wraps_origin: bool = False


@dataclass(frozen=True, slots=True)
class EnzymeSiteOverlay:
    """Map-only restriction site overlay that can represent circular wrapping sites."""

    enzyme_name: str
    recognition_sequence: str
    start: int
    end: int
    strand: int = 1
    wraps_origin: bool = False


@dataclass(frozen=True, slots=True)
class PrimerBindingOverlay:
    """Map-only primer binding overlay that can represent circular wrapping bindings."""

    name: str
    sequence: str
    start: int
    end: int
    strand: int
    wraps_origin: bool = False


@dataclass(frozen=True, slots=True)
class MapOverlays:
    """Biological overlays prepared outside the GUI for map rendering."""

    enzyme_sites: tuple[EnzymeSiteOverlay, ...] = ()
    primer_bindings: tuple[PrimerBindingOverlay, ...] = ()


@dataclass(frozen=True, slots=True)
class Label:
    """A positioned text label."""

    text: str
    position: Point
    anchor: str = "middle"
    color: str = "#212529"


@dataclass(frozen=True, slots=True)
class ScaleBar:
    """Map scale annotation."""

    start_point: Point
    end_point: Point
    length_bp: int
    label: str


@dataclass(frozen=True, slots=True)
class OriginMarker:
    """Sequence origin marker."""

    position: int
    start_point: Point
    end_point: Point
    label: str = "0"
    angle: float | None = None


@dataclass(frozen=True, slots=True)
class PlasmidMapModel:
    """Intermediate drawing model for plasmid maps."""

    record_id: str
    sequence_length: int
    topology: MoleculeTopology
    width: int
    height: int
    map_kind: MapKind
    center: Point | None
    radius: float | None
    backbone_start: Point | None
    backbone_end: Point | None
    feature_arcs: tuple[FeatureArc, ...]
    feature_arrows: tuple[FeatureArrow, ...]
    enzyme_ticks: tuple[EnzymeTick, ...]
    primer_arrows: tuple[PrimerArrow, ...]
    labels: tuple[Label, ...]
    scale: ScaleBar
    origin_marker: OriginMarker


def render_plasmid_map(
    record: SequenceRecord,
    *,
    overlays: MapOverlays | None = None,
    style: MapStyle | None = None,
    feature_colors: Mapping[str, str] | None = None,
) -> PlasmidMapModel:
    """Create a GUI-independent drawing model for a sequence record."""

    if record.length == 0:
        msg = "cannot render an empty sequence"
        raise ValueError(msg)
    style = style or MapStyle()
    colors = dict(style.feature_colors)
    if feature_colors:
        colors.update(feature_colors)

    overlay_model = overlays or prepare_map_overlays(record)

    if record.topology is MoleculeTopology.CIRCULAR:
        return _render_circular(record, overlay_model, style, colors)
    return _render_linear(record, overlay_model, style, colors)


def prepare_map_overlays(
    record: SequenceRecord,
    *,
    restriction_analysis: RestrictionAnalysis | None = None,
    primer_bindings_by_primer_index: Mapping[int, Iterable[PrimerBinding]] | None = None,
) -> MapOverlays:
    """Prepare map overlays from core analysis results without GUI coordinate logic."""

    enzyme_overlays = [
        _enzyme_overlay_from_enzyme_site(site)
        for site in record.enzyme_sites
    ]
    if restriction_analysis is not None:
        enzyme_overlays.extend(
            _enzyme_overlay_from_restriction_site(site)
            for site in restriction_analysis.sites
        )

    primer_overlays = [
        _primer_overlay_from_primer(primer)
        for primer in record.primers
        if primer.start is not None and primer.end is not None
    ]
    if primer_bindings_by_primer_index:
        for primer_index, bindings in primer_bindings_by_primer_index.items():
            if primer_index < 0 or primer_index >= len(record.primers):
                continue
            primer = record.primers[primer_index]
            primer_overlays.extend(
                _primer_overlay_from_binding(primer.name, primer.sequence, binding)
                for binding in bindings
            )

    return MapOverlays(
        enzyme_sites=tuple(enzyme_overlays),
        primer_bindings=tuple(primer_overlays),
    )


def _enzyme_overlay_from_enzyme_site(site: EnzymeSite) -> EnzymeSiteOverlay:
    return EnzymeSiteOverlay(
        enzyme_name=site.enzyme_name,
        recognition_sequence=site.recognition_sequence,
        start=site.start,
        end=site.end,
        strand=site.strand,
        wraps_origin=False,
    )


def _enzyme_overlay_from_restriction_site(site: RestrictionSite) -> EnzymeSiteOverlay:
    return EnzymeSiteOverlay(
        enzyme_name=site.enzyme_name,
        recognition_sequence=site.recognition_sequence,
        start=site.start,
        end=site.end,
        strand=site.strand,
        wraps_origin=site.wraps_origin,
    )


def _primer_overlay_from_primer(primer: Primer) -> PrimerBindingOverlay:
    assert primer.start is not None
    assert primer.end is not None
    return PrimerBindingOverlay(
        name=primer.name,
        sequence=primer.sequence,
        start=primer.start,
        end=primer.end,
        strand=primer.strand,
        wraps_origin=False,
    )


def _primer_overlay_from_binding(
    primer_name: str,
    primer_sequence: str,
    binding: PrimerBinding,
) -> PrimerBindingOverlay:
    return PrimerBindingOverlay(
        name=primer_name,
        sequence=primer_sequence,
        start=binding.start,
        end=binding.end,
        strand=binding.strand,
        wraps_origin=binding.wraps_origin,
    )


def to_svg(model: PlasmidMapModel, *, background: str | None = None) -> str:
    """Export a plasmid map drawing model as SVG text."""

    background_color = background or "#ffffff"
    parts: list[str] = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{model.width}" '
            f'height="{model.height}" viewBox="0 0 {model.width} {model.height}" '
            'role="img">'
        ),
        f"<title>{escape(model.record_id)} plasmid map</title>",
        f'<rect width="100%" height="100%" fill="{background_color}"/>',
    ]
    if model.map_kind == "circular" and model.center and model.radius is not None:
        parts.append(
            f'<circle cx="{_fmt(model.center.x)}" cy="{_fmt(model.center.y)}" '
            f'r="{_fmt(model.radius)}" fill="none" stroke="#343a40" stroke-width="2"/>'
        )
    elif model.backbone_start and model.backbone_end:
        parts.append(
            f'<line x1="{_fmt(model.backbone_start.x)}" y1="{_fmt(model.backbone_start.y)}" '
            f'x2="{_fmt(model.backbone_end.x)}" y2="{_fmt(model.backbone_end.y)}" '
            'stroke="#343a40" stroke-width="3" stroke-linecap="round"/>'
        )

    for arc in model.feature_arcs:
        if arc.map_kind == "circular" and arc.center and arc.radius is not None:
            parts.append(_svg_arc(arc))
        elif arc.start_point and arc.end_point:
            parts.append(
                f'<line x1="{_fmt(arc.start_point.x)}" y1="{_fmt(arc.start_point.y)}" '
                f'x2="{_fmt(arc.end_point.x)}" y2="{_fmt(arc.end_point.y)}" '
                f'stroke="{arc.color}" stroke-width="12" stroke-linecap="round" fill="none"/>'
            )

    for arrow in model.feature_arrows:
        parts.append(_svg_polygon(arrow.points, arrow.color))
    for primer in model.primer_arrows:
        parts.append(_svg_polygon(primer.points, primer.color))
    for tick in model.enzyme_ticks:
        parts.append(
            f'<line x1="{_fmt(tick.start_point.x)}" y1="{_fmt(tick.start_point.y)}" '
            f'x2="{_fmt(tick.end_point.x)}" y2="{_fmt(tick.end_point.y)}" '
            'stroke="#212529" stroke-width="2"/>'
        )

    origin = model.origin_marker
    parts.append(
        f'<line x1="{_fmt(origin.start_point.x)}" y1="{_fmt(origin.start_point.y)}" '
        f'x2="{_fmt(origin.end_point.x)}" y2="{_fmt(origin.end_point.y)}" '
        'stroke="#d9480f" stroke-width="3"/>'
    )
    parts.append(
        f'<line x1="{_fmt(model.scale.start_point.x)}" y1="{_fmt(model.scale.start_point.y)}" '
        f'x2="{_fmt(model.scale.end_point.x)}" y2="{_fmt(model.scale.end_point.y)}" '
        'stroke="#495057" stroke-width="2"/>'
    )

    for label in model.labels:
        parts.append(
            f'<text x="{_fmt(label.position.x)}" y="{_fmt(label.position.y)}" '
            f'text-anchor="{escape(label.anchor)}" font-family="Arial, sans-serif" '
            f'font-size="12" fill="{label.color}">{escape(label.text)}</text>'
        )
    parts.append("</svg>")
    return "\n".join(parts)


def write_svg(model: PlasmidMapModel, path: str | Path) -> None:
    """Write SVG output to a path."""

    Path(path).write_text(to_svg(model), encoding="utf-8")


def _render_circular(
    record: SequenceRecord,
    overlays: MapOverlays,
    style: MapStyle,
    colors: Mapping[str, str],
) -> PlasmidMapModel:
    center = Point(style.width / 2, style.height / 2)
    radius = style.circular_radius
    lanes = _assign_feature_lanes(record.features)
    feature_arcs: list[FeatureArc] = []
    feature_arrows: list[FeatureArrow] = []
    labels: list[Label] = []

    for feature_index, feature in enumerate(record.features):
        lane = lanes[feature_index]
        feature_radius = radius + style.feature_lane_gap * (lane + 1)
        color = _feature_color(feature, colors)
        label = _feature_label(feature)
        for segment_index, segment in enumerate(feature.segments):
            start_angle = _angle(segment.start, record.length)
            end_angle = _angle(segment.end, record.length)
            start_point = _polar(center, feature_radius, start_angle)
            end_point = _polar(center, feature_radius, end_angle)
            feature_arcs.append(
                FeatureArc(
                    feature_index=feature_index,
                    segment_index=segment_index,
                    feature_type=feature.type,
                    label=label,
                    start=segment.start,
                    end=segment.end,
                    strand=segment.strand,
                    lane=lane,
                    color=color,
                    map_kind="circular",
                    start_angle=start_angle,
                    end_angle=end_angle,
                    radius=feature_radius,
                    center=center,
                    start_point=start_point,
                    end_point=end_point,
                )
            )
            if segment.strand:
                feature_arrows.append(
                    _circular_arrow(
                        feature_index,
                        segment_index,
                        segment,
                        center,
                        feature_radius,
                        record.length,
                        color,
                    )
                )
            labels.append(
                Label(
                    text=label,
                    position=_polar(
                        center,
                        feature_radius + 22,
                        _angle(_segment_midpoint(segment, record.length), record.length),
                    ),
                )
            )

    enzyme_ticks = tuple(_circular_enzyme_ticks(overlays.enzyme_sites, center, radius, record.length))
    primer_arrows = tuple(
        _circular_primer_arrow(
            primer,
            center,
            radius + style.feature_lane_gap * (max(lanes, default=-1) + 3),
            record.length,
        )
        for primer in overlays.primer_bindings
    )
    labels.extend(_enzyme_labels(enzyme_ticks, center, radius + 42))
    labels.extend(_primer_labels(primer_arrows))

    origin_marker = OriginMarker(
        position=0,
        start_point=_polar(center, radius - 15, -90),
        end_point=_polar(center, radius + 20, -90),
        angle=-90,
    )
    labels.append(Label(text="0", position=_polar(center, radius + 35, -90)))
    scale = _scale_bar(style, record.length)
    labels.append(Label(text=scale.label, position=Point((scale.start_point.x + scale.end_point.x) / 2, scale.start_point.y - 8)))

    return PlasmidMapModel(
        record_id=record.id,
        sequence_length=record.length,
        topology=record.topology,
        width=style.width,
        height=style.height,
        map_kind="circular",
        center=center,
        radius=radius,
        backbone_start=None,
        backbone_end=None,
        feature_arcs=tuple(feature_arcs),
        feature_arrows=tuple(feature_arrows),
        enzyme_ticks=enzyme_ticks,
        primer_arrows=primer_arrows,
        labels=tuple(labels),
        scale=scale,
        origin_marker=origin_marker,
    )


def _render_linear(
    record: SequenceRecord,
    overlays: MapOverlays,
    style: MapStyle,
    colors: Mapping[str, str],
) -> PlasmidMapModel:
    y = style.height / 2
    backbone_start = Point(style.margin, y)
    backbone_end = Point(style.width - style.margin, y)
    lanes = _assign_feature_lanes(record.features)
    feature_arcs: list[FeatureArc] = []
    feature_arrows: list[FeatureArrow] = []
    labels: list[Label] = []

    for feature_index, feature in enumerate(record.features):
        lane = lanes[feature_index]
        feature_y = y - 35 - lane * style.feature_lane_gap
        color = _feature_color(feature, colors)
        label = _feature_label(feature)
        for segment_index, segment in enumerate(feature.segments):
            start_point = Point(_linear_x(segment.start, record.length, style), feature_y)
            end_point = Point(_linear_x(segment.end, record.length, style), feature_y)
            feature_arcs.append(
                FeatureArc(
                    feature_index=feature_index,
                    segment_index=segment_index,
                    feature_type=feature.type,
                    label=label,
                    start=segment.start,
                    end=segment.end,
                    strand=segment.strand,
                    lane=lane,
                    color=color,
                    map_kind="linear",
                    start_point=start_point,
                    end_point=end_point,
                )
            )
            if segment.strand:
                feature_arrows.append(_linear_feature_arrow(feature_index, segment_index, segment, feature_y, record.length, style, color))
            labels.append(Label(text=label, position=Point((start_point.x + end_point.x) / 2, feature_y - 16)))

    enzyme_ticks = tuple(_linear_enzyme_ticks(overlays.enzyme_sites, style, y, record.length))
    primer_y = y + 44
    primer_arrows = tuple(
        _linear_primer_arrow(primer, primer_y, record.length, style)
        for primer in overlays.primer_bindings
    )
    labels.extend(_linear_enzyme_labels(enzyme_ticks))
    labels.extend(_primer_labels(primer_arrows))

    origin_marker = OriginMarker(
        position=0,
        start_point=Point(style.margin, y - 18),
        end_point=Point(style.margin, y + 18),
    )
    labels.append(Label(text="0", position=Point(style.margin, y + 36)))
    scale = _scale_bar(style, record.length)
    labels.append(Label(text=scale.label, position=Point((scale.start_point.x + scale.end_point.x) / 2, scale.start_point.y - 8)))

    return PlasmidMapModel(
        record_id=record.id,
        sequence_length=record.length,
        topology=record.topology,
        width=style.width,
        height=style.height,
        map_kind="linear",
        center=None,
        radius=None,
        backbone_start=backbone_start,
        backbone_end=backbone_end,
        feature_arcs=tuple(feature_arcs),
        feature_arrows=tuple(feature_arrows),
        enzyme_ticks=enzyme_ticks,
        primer_arrows=primer_arrows,
        labels=tuple(labels),
        scale=scale,
        origin_marker=origin_marker,
    )


def _assign_feature_lanes(features: Iterable[Feature]) -> tuple[int, ...]:
    lane_intervals: list[list[tuple[int, int]]] = []
    assignments: list[int] = []
    for feature in features:
        intervals = [(segment.start, segment.end) for segment in feature.segments]
        for lane, occupied in enumerate(lane_intervals):
            if not any(_intervals_overlap(interval, used) for interval in intervals for used in occupied):
                occupied.extend(intervals)
                assignments.append(lane)
                break
        else:
            lane_intervals.append(list(intervals))
            assignments.append(len(lane_intervals) - 1)
    return tuple(assignments)


def _intervals_overlap(left: tuple[int, int], right: tuple[int, int]) -> bool:
    return left[0] < right[1] and right[0] < left[1]


def _feature_color(feature: Feature, colors: Mapping[str, str]) -> str:
    return colors.get(feature.type, colors.get("misc_feature", "#868e96"))


def _feature_label(feature: Feature) -> str:
    for key in ("label", "gene", "product", "note"):
        value = feature.qualifiers.get(key)
        if value:
            return _qualifier_text(value)
    return feature.name or feature.type


def _qualifier_text(value: object) -> str:
    if isinstance(value, tuple | list):
        for item in value:
            text = str(item).strip()
            if text:
                return text
        return ""
    return str(value)


def _angle(position: float, sequence_length: int) -> float:
    return -90.0 + (position / sequence_length) * 360.0


def _polar(center: Point, radius: float, angle_degrees: float) -> Point:
    radians = angle_degrees * pi / 180.0
    return Point(center.x + radius * cos(radians), center.y + radius * sin(radians))


def _segment_midpoint(segment: FeatureSegment, sequence_length: int) -> float:
    return (segment.start + segment.end) / 2 % sequence_length


def _linear_x(position: int, sequence_length: int, style: MapStyle) -> float:
    usable = style.width - 2 * style.margin
    return style.margin + usable * (position / sequence_length)


def _circular_arrow(
    feature_index: int,
    segment_index: int,
    segment: FeatureSegment,
    center: Point,
    radius: float,
    sequence_length: int,
    color: str,
) -> FeatureArrow:
    tip_position = segment.end if segment.strand == 1 else segment.start
    angle = _angle(tip_position, sequence_length)
    direction = 1 if segment.strand == 1 else -1
    tip = _polar(center, radius, angle)
    wing_a = _polar(center, radius - 8, angle - direction * 5)
    wing_b = _polar(center, radius + 8, angle - direction * 5)
    return FeatureArrow(
        feature_index=feature_index,
        segment_index=segment_index,
        strand=segment.strand,
        color=color,
        points=(tip, wing_a, wing_b),
    )


def _linear_feature_arrow(
    feature_index: int,
    segment_index: int,
    segment: FeatureSegment,
    y: float,
    sequence_length: int,
    style: MapStyle,
    color: str,
) -> FeatureArrow:
    tip_x = _linear_x(segment.end if segment.strand == 1 else segment.start, sequence_length, style)
    direction = 1 if segment.strand == 1 else -1
    return FeatureArrow(
        feature_index=feature_index,
        segment_index=segment_index,
        strand=segment.strand,
        color=color,
        points=(
            Point(tip_x, y),
            Point(tip_x - direction * 14, y - 8),
            Point(tip_x - direction * 14, y + 8),
        ),
    )


def _circular_enzyme_ticks(
    enzyme_sites: Iterable[EnzymeSiteOverlay],
    center: Point,
    radius: float,
    sequence_length: int,
) -> Iterable[EnzymeTick]:
    for site in enzyme_sites:
        angle = _angle(site.start, sequence_length)
        yield EnzymeTick(
            enzyme_name=site.enzyme_name,
            position=site.start,
            recognition_sequence=site.recognition_sequence,
            start_point=_polar(center, radius - 10, angle),
            end_point=_polar(center, radius + 18, angle),
            angle=angle,
            site_start=site.start,
            site_end=site.end,
            wraps_origin=site.wraps_origin,
        )


def _linear_enzyme_ticks(
    enzyme_sites: Iterable[EnzymeSiteOverlay],
    style: MapStyle,
    y: float,
    sequence_length: int,
) -> Iterable[EnzymeTick]:
    for site in enzyme_sites:
        if site.wraps_origin or site.end <= site.start:
            continue
        x = _linear_x(site.start, sequence_length, style)
        yield EnzymeTick(
            enzyme_name=site.enzyme_name,
            position=site.start,
            recognition_sequence=site.recognition_sequence,
            start_point=Point(x, y - 15),
            end_point=Point(x, y + 15),
            site_start=site.start,
            site_end=site.end,
            wraps_origin=site.wraps_origin,
        )


def _circular_primer_arrow(
    primer: PrimerBindingOverlay,
    center: Point,
    radius: float,
    sequence_length: int,
) -> PrimerArrow:
    tip_position = primer.end if primer.strand == 1 else primer.start
    angle = _angle(tip_position, sequence_length)
    direction = 1 if primer.strand == 1 else -1
    tip = _polar(center, radius, angle)
    return PrimerArrow(
        name=primer.name,
        start=primer.start,
        end=primer.end,
        strand=primer.strand,
        color=DEFAULT_FEATURE_COLORS["primer"],
        points=(tip, _polar(center, radius - 7, angle - direction * 5), _polar(center, radius + 7, angle - direction * 5)),
        wraps_origin=primer.wraps_origin,
    )


def _linear_primer_arrow(
    primer: PrimerBindingOverlay,
    y: float,
    sequence_length: int,
    style: MapStyle,
) -> PrimerArrow:
    if primer.end <= primer.start:
        msg = "linear primer overlays must not wrap"
        raise ValueError(msg)
    tip_x = _linear_x(primer.end if primer.strand == 1 else primer.start, sequence_length, style)
    direction = 1 if primer.strand == 1 else -1
    return PrimerArrow(
        name=primer.name,
        start=primer.start,
        end=primer.end,
        strand=primer.strand,
        color=DEFAULT_FEATURE_COLORS["primer"],
        points=(
            Point(tip_x, y),
            Point(tip_x - direction * 14, y - 8),
            Point(tip_x - direction * 14, y + 8),
        ),
        wraps_origin=primer.wraps_origin,
    )


def _enzyme_labels(ticks: Iterable[EnzymeTick], center: Point, radius: float) -> list[Label]:
    labels: list[Label] = []
    for tick in ticks:
        angle = tick.angle if tick.angle is not None else -90
        labels.append(Label(text=tick.enzyme_name, position=_polar(center, radius, angle)))
    return labels


def _linear_enzyme_labels(ticks: Iterable[EnzymeTick]) -> list[Label]:
    return [Label(text=tick.enzyme_name, position=Point(tick.end_point.x, tick.end_point.y + 18)) for tick in ticks]


def _primer_labels(primers: Iterable[PrimerArrow]) -> list[Label]:
    labels: list[Label] = []
    for primer in primers:
        tip = primer.points[0]
        labels.append(Label(text=primer.name, position=Point(tip.x, tip.y + 18)))
    return labels


def _scale_bar(style: MapStyle, sequence_length: int) -> ScaleBar:
    length_bp = _nice_scale_length(sequence_length)
    pixel_length = min(180, (style.width - 2 * style.margin) * length_bp / sequence_length)
    start = Point(style.margin, style.height - style.margin)
    end = Point(style.margin + pixel_length, style.height - style.margin)
    return ScaleBar(start_point=start, end_point=end, length_bp=length_bp, label=f"{length_bp} bp")


def _nice_scale_length(sequence_length: int) -> int:
    candidates = (10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000)
    target = max(1, sequence_length // 5)
    return max(candidate for candidate in candidates if candidate <= max(target, candidates[0]))


def _svg_arc(arc: FeatureArc) -> str:
    assert arc.center is not None
    assert arc.radius is not None
    assert arc.start_point is not None
    assert arc.end_point is not None
    assert arc.start_angle is not None
    assert arc.end_angle is not None
    angle_delta = (arc.end_angle - arc.start_angle) % 360
    if angle_delta == 0 and arc.end > arc.start:
        return (
            f'<circle cx="{_fmt(arc.center.x)}" cy="{_fmt(arc.center.y)}" '
            f'r="{_fmt(arc.radius)}" fill="none" stroke="{arc.color}" '
            'stroke-width="12"/>'
        )
    large_arc = 1 if angle_delta > 180 else 0
    return (
        f'<path d="M {_fmt(arc.start_point.x)} {_fmt(arc.start_point.y)} '
        f'A {_fmt(arc.radius)} {_fmt(arc.radius)} 0 {large_arc} 1 '
        f'{_fmt(arc.end_point.x)} {_fmt(arc.end_point.y)}" '
        f'stroke="{arc.color}" stroke-width="12" stroke-linecap="round" fill="none"/>'
    )


def _svg_polygon(points: tuple[Point, Point, Point], color: str) -> str:
    point_text = " ".join(f"{_fmt(point.x)},{_fmt(point.y)}" for point in points)
    return f'<polygon points="{point_text}" fill="{color}"/>'


def _fmt(value: float) -> str:
    return f"{value:.3f}".rstrip("0").rstrip(".")
