"""PySide6 widget for approximate agarose gel simulation models."""

from __future__ import annotations

import sys
from pathlib import Path

from plasmidlab.core.gel import AgaroseGelModel, GelLaneInput, simulate_gel
from plasmidlab.render.agarose_gel import GelRenderStyle, gel_to_svg, write_gel_svg

try:
    from PySide6.QtCore import QPointF, QRectF, Qt
    from PySide6.QtGui import QColor, QFont, QImage, QPainter, QPen
    from PySide6.QtWidgets import QApplication, QWidget
except ImportError:  # pragma: no cover - exercised only without PySide6.
    QApplication = None  # type: ignore[assignment]
    QWidget = object  # type: ignore[assignment,misc]
    QPointF = None  # type: ignore[assignment]
    QRectF = None  # type: ignore[assignment]
    Qt = None  # type: ignore[assignment]
    QColor = None  # type: ignore[assignment]
    QFont = None  # type: ignore[assignment]
    QImage = None  # type: ignore[assignment]
    QPainter = None  # type: ignore[assignment]
    QPen = None  # type: ignore[assignment]


class AgaroseGelWidget(QWidget):
    """Render an approximate agarose gel model using PySide6."""

    def __init__(
        self,
        model: AgaroseGelModel | None = None,
        *,
        style: GelRenderStyle | None = None,
        parent: object | None = None,
    ) -> None:
        if QApplication is None:
            msg = "PySide6 is required to use AgaroseGelWidget"
            raise ImportError(msg)
        super().__init__(parent)  # type: ignore[misc]
        self._model = model or sample_gel_model()
        self._style = style or GelRenderStyle()
        self.setMinimumSize(self._style.width, self._style.height)
        self.setWindowTitle("Agarose gel")

    @property
    def model(self) -> AgaroseGelModel:
        """Return the gel model currently rendered by the widget."""

        return self._model

    def set_model(self, model: AgaroseGelModel) -> None:
        """Replace the gel model and schedule repaint."""

        self._model = model
        self.update()

    def to_svg(self) -> str:
        """Return SVG output for the current model."""

        return gel_to_svg(self._model, style=self._style)

    def write_svg(self, path: str | Path) -> None:
        """Write the current gel model as SVG."""

        write_gel_svg(self._model, path, style=self._style)

    def write_png(self, path: str | Path) -> None:
        """Render the current widget to a PNG image."""

        image = QImage(self._style.width, self._style.height, QImage.Format.Format_ARGB32)
        image.fill(QColor("#ffffff"))
        self.resize(self._style.width, self._style.height)
        self.render(image)
        image.save(str(path), "PNG")

    def paintEvent(self, event: object) -> None:  # noqa: N802 - Qt API name.
        del event
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.fillRect(self.rect(), QColor(self._style.background_color))
        layout = _layout(self._model, self._style)
        gel_left, gel_top, gel_width, gel_height = layout["gel_rect"]
        painter.setPen(QPen(QColor(self._style.gel_border_color), 1))
        painter.setBrush(QColor(self._style.gel_color))
        painter.drawRoundedRect(QRectF(gel_left, gel_top, gel_width, gel_height), 4, 4)
        self._draw_lanes(painter, layout)
        self._draw_footer(painter)
        painter.end()

    def _draw_lanes(self, painter: object, layout: dict[str, object]) -> None:
        painter.setFont(QFont("Arial", 10))
        for lane_index, lane in enumerate(self._model.lanes):
            lane_center = layout["lane_centers"][lane_index]
            well_x = lane_center - self._style.lane_width / 2
            painter.setPen(QPen(QColor("#91a7ff"), 1))
            painter.setBrush(QColor("#ffffff"))
            painter.drawRoundedRect(
                QRectF(well_x, layout["well_y"], self._style.lane_width, self._style.well_height),
                2,
                2,
            )
            painter.setPen(QPen(QColor(self._style.label_color), 1))
            painter.drawText(
                QRectF(lane_center - 80, self._style.margin_top - 44, 160, 24),
                Qt.AlignmentFlag.AlignCenter,
                lane.name,
            )
            for band in lane.bands:
                y = float(layout["run_start_y"]) + band.migration * float(layout["run_height"])
                intensity = max(0.0, min(1.0, band.relative_intensity))
                width = self._style.lane_width * (0.42 + 0.58 * intensity)
                x = lane_center - width / 2
                alpha = int(75 + 180 * intensity)
                band_color = QColor(self._style.ladder_band_color if lane.is_ladder else self._style.band_color)
                band_color.setAlpha(alpha)
                painter.setPen(Qt.PenStyle.NoPen)
                painter.setBrush(band_color)
                painter.drawRoundedRect(QRectF(x, y, width, self._style.band_height), 2, 2)
                painter.setPen(QPen(QColor(self._style.label_color), 1))
                painter.drawText(QPointF(lane_center + self._style.lane_width / 2 + 8, y + 9), band.label)

    def _draw_footer(self, painter: object) -> None:
        painter.setFont(QFont("Arial", 9))
        painter.setPen(QPen(QColor("#495057"), 1))
        painter.drawText(
            QRectF(0, self._style.height - 42, self._style.width, 24),
            Qt.AlignmentFlag.AlignCenter,
            (
                f"{self._model.agarose_percentage:g}% agarose, {self._model.voltage:g} V, "
                f"{self._model.run_time_minutes:g} min, {self._model.amount_mode.value}"
            ),
        )


def sample_gel_model() -> AgaroseGelModel:
    """Return a deterministic synthetic gel for demos and smoke tests."""

    return simulate_gel(
        (
            GelLaneInput(name="Digest", fragments=(3000, 1500, 750)),
            GelLaneInput(name="PCR", fragments=(850,)),
        ),
        ladder="1kb_dna_ladder",
        amount_mode="equal_mass",
    )


def show_sample_gel() -> int:
    """Open a standalone sample gel window."""

    if QApplication is None:
        msg = "PySide6 is required to open the sample gel widget"
        raise ImportError(msg)
    app = QApplication.instance() or QApplication(sys.argv)
    widget = AgaroseGelWidget(sample_gel_model())
    widget.show()
    return app.exec()


def _layout(model: AgaroseGelModel, style: GelRenderStyle) -> dict[str, object]:
    lane_count = max(1, len(model.lanes))
    gel_left = style.margin_left
    gel_top = style.margin_top
    gel_width = style.width - style.margin_left - style.margin_right
    gel_height = style.height - style.margin_top - style.margin_bottom
    lane_spacing = gel_width / lane_count
    return {
        "gel_rect": (gel_left, gel_top, gel_width, gel_height),
        "well_y": gel_top + 16,
        "run_start_y": gel_top + 38,
        "run_height": gel_height - 58,
        "lane_centers": tuple(gel_left + lane_spacing * (index + 0.5) for index in range(lane_count)),
    }


if __name__ == "__main__":
    raise SystemExit(show_sample_gel())
