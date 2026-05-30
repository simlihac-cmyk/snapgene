"""GUI-independent agarose gel SVG rendering."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from xml.sax.saxutils import escape

from plasmidlab.core.gel import AgaroseGelModel, GelBand


@dataclass(frozen=True, slots=True)
class GelRenderStyle:
    """Layout settings for agarose gel rendering."""

    width: int = 900
    height: int = 620
    margin_left: int = 90
    margin_right: int = 130
    margin_top: int = 70
    margin_bottom: int = 70
    lane_width: int = 72
    well_height: int = 16
    band_height: int = 7
    background_color: str = "#ffffff"
    gel_color: str = "#e7f5ff"
    gel_border_color: str = "#74c0fc"
    band_color: str = "#212529"
    ladder_band_color: str = "#1c7ed6"
    label_color: str = "#212529"


def gel_to_svg(model: AgaroseGelModel, *, style: GelRenderStyle | None = None) -> str:
    """Export an agarose gel model as SVG text."""

    style = style or GelRenderStyle()
    layout = _layout(model, style)
    gel_left, gel_top, gel_width, gel_height = layout["gel_rect"]
    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{style.width}" '
            f'height="{style.height}" viewBox="0 0 {style.width} {style.height}" role="img">'
        ),
        "<title>Agarose gel simulation</title>",
        f'<rect width="100%" height="100%" fill="{style.background_color}"/>',
        (
            f'<rect x="{_fmt(gel_left)}" y="{_fmt(gel_top)}" width="{_fmt(gel_width)}" '
            f'height="{_fmt(gel_height)}" rx="4" fill="{style.gel_color}" '
            f'stroke="{style.gel_border_color}" stroke-width="1"/>'
        ),
    ]
    for lane_index, lane in enumerate(model.lanes):
        lane_center = layout["lane_centers"][lane_index]
        well_x = lane_center - style.lane_width / 2
        parts.append(
            f'<rect x="{_fmt(well_x)}" y="{_fmt(layout["well_y"])}" '
            f'width="{style.lane_width}" height="{style.well_height}" rx="2" '
            'fill="#ffffff" stroke="#91a7ff" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{_fmt(lane_center)}" y="{_fmt(style.margin_top - 22)}" '
            f'text-anchor="middle" font-family="Arial, sans-serif" font-size="12" '
            f'fill="{style.label_color}">{escape(lane.name)}</text>'
        )
        for band in lane.bands:
            parts.append(_svg_band(band, lane_center, lane.is_ladder, layout, style))
    parts.append(
        f'<text x="{_fmt(style.width / 2)}" y="{_fmt(style.height - 24)}" text-anchor="middle" '
        f'font-family="Arial, sans-serif" font-size="11" fill="#495057">'
        f'{model.agarose_percentage:g}% agarose, {model.voltage:g} V, '
        f'{model.run_time_minutes:g} min, {escape(model.amount_mode.value)}</text>'
    )
    parts.append("</svg>")
    return "\n".join(parts)


def write_gel_svg(
    model: AgaroseGelModel,
    path: str | Path,
    *,
    style: GelRenderStyle | None = None,
) -> None:
    """Write agarose gel SVG output to a path."""

    Path(path).write_text(gel_to_svg(model, style=style), encoding="utf-8")


def _svg_band(
    band: GelBand,
    lane_center: float,
    is_ladder: bool,
    layout: dict[str, object],
    style: GelRenderStyle,
) -> str:
    y = _band_y(band, layout)
    width = style.lane_width * (0.42 + 0.58 * max(0.0, min(1.0, band.relative_intensity)))
    x = lane_center - width / 2
    opacity = 0.30 + 0.70 * max(0.0, min(1.0, band.relative_intensity))
    color = style.ladder_band_color if is_ladder else style.band_color
    label_x = lane_center + style.lane_width / 2 + 8
    return "\n".join(
        (
            (
                f'<rect x="{_fmt(x)}" y="{_fmt(y)}" width="{_fmt(width)}" '
                f'height="{style.band_height}" rx="2" fill="{color}" opacity="{_fmt(opacity)}"/>'
            ),
            (
                f'<text x="{_fmt(label_x)}" y="{_fmt(y + style.band_height + 3)}" '
                f'text-anchor="start" font-family="Arial, sans-serif" font-size="10" '
                f'fill="{style.label_color}">{escape(band.label)}</text>'
            ),
        )
    )


def _layout(model: AgaroseGelModel, style: GelRenderStyle) -> dict[str, object]:
    lane_count = max(1, len(model.lanes))
    gel_left = style.margin_left
    gel_top = style.margin_top
    gel_width = style.width - style.margin_left - style.margin_right
    gel_height = style.height - style.margin_top - style.margin_bottom
    lane_spacing = gel_width / lane_count
    lane_centers = tuple(gel_left + lane_spacing * (index + 0.5) for index in range(lane_count))
    return {
        "gel_rect": (gel_left, gel_top, gel_width, gel_height),
        "well_y": gel_top + 16,
        "run_start_y": gel_top + 38,
        "run_height": gel_height - 58,
        "lane_centers": lane_centers,
    }


def _band_y(band: GelBand, layout: dict[str, object]) -> float:
    return float(layout["run_start_y"]) + band.migration * float(layout["run_height"])


def _fmt(value: float) -> str:
    return f"{value:.3f}".rstrip("0").rstrip(".")
