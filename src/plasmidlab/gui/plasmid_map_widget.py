"""Minimal PySide6 widget for rendering a PlasmidLab plasmid map model."""

from __future__ import annotations

import sys

from plasmidlab.core import EnzymeSite, Feature, MoleculeTopology, Primer, SequenceRecord
from plasmidlab.render import FeatureArc, MapStyle, PlasmidMapModel, Point, render_plasmid_map

try:
    from PySide6.QtCore import QPointF, QRectF, Qt
    from PySide6.QtGui import QBrush, QColor, QFont, QPainter, QPen, QPolygonF
    from PySide6.QtWidgets import QApplication, QWidget
except ImportError:  # pragma: no cover - exercised only in environments without PySide6.
    QApplication = None  # type: ignore[assignment]
    QWidget = object  # type: ignore[assignment,misc]
    QPainter = None  # type: ignore[assignment]
    QPen = None  # type: ignore[assignment]
    QBrush = None  # type: ignore[assignment]
    QColor = None  # type: ignore[assignment]
    QPointF = None  # type: ignore[assignment]
    QRectF = None  # type: ignore[assignment]
    Qt = None  # type: ignore[assignment]
    QFont = None  # type: ignore[assignment]
    QPolygonF = None  # type: ignore[assignment]


class PlasmidMapWidget(QWidget):
    """Render a PlasmidMapModel using PySide6 painting primitives."""

    def __init__(
        self,
        record: SequenceRecord | None = None,
        *,
        model: PlasmidMapModel | None = None,
        style: MapStyle | None = None,
        parent: object | None = None,
    ) -> None:
        if QApplication is None:
            msg = "PySide6 is required to use PlasmidMapWidget"
            raise ImportError(msg)
        super().__init__(parent)  # type: ignore[misc]
        if model is None:
            model = render_plasmid_map(record or sample_plasmid_record(), style=style)
        self._model = model
        self.setMinimumSize(model.width, model.height)
        self.setWindowTitle(f"{model.record_id} plasmid map")

    @property
    def model(self) -> PlasmidMapModel:
        """Return the drawing model currently rendered by the widget."""

        return self._model

    def set_model(self, model: PlasmidMapModel) -> None:
        """Replace the drawing model and schedule repaint."""

        self._model = model
        self.setMinimumSize(model.width, model.height)
        self.update()

    def paintEvent(self, event: object) -> None:  # noqa: N802 - Qt API name.
        del event
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.fillRect(self.rect(), QColor("#ffffff"))
        self._draw_backbone(painter)
        self._draw_features(painter)
        self._draw_ticks(painter)
        self._draw_arrows(painter)
        self._draw_scale_and_origin(painter)
        self._draw_labels(painter)
        painter.end()

    def _draw_backbone(self, painter: object) -> None:
        pen = QPen(QColor("#343a40"), 2)
        painter.setPen(pen)
        painter.setBrush(Qt.BrushStyle.NoBrush)
        model = self._model
        if model.map_kind == "circular" and model.center and model.radius is not None:
            painter.drawEllipse(
                QPointF(model.center.x, model.center.y),
                model.radius,
                model.radius,
            )
        elif model.backbone_start and model.backbone_end:
            painter.drawLine(_qpoint(model.backbone_start), _qpoint(model.backbone_end))

    def _draw_features(self, painter: object) -> None:
        for arc in self._model.feature_arcs:
            pen = QPen(QColor(arc.color), 12)
            pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            painter.setPen(pen)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            if arc.map_kind == "circular":
                self._draw_circular_arc(painter, arc)
            elif arc.start_point and arc.end_point:
                painter.drawLine(_qpoint(arc.start_point), _qpoint(arc.end_point))

    def _draw_circular_arc(self, painter: object, arc: FeatureArc) -> None:
        assert arc.center is not None
        assert arc.radius is not None
        assert arc.start_angle is not None
        assert arc.end_angle is not None
        rect = QRectF(
            arc.center.x - arc.radius,
            arc.center.y - arc.radius,
            arc.radius * 2,
            arc.radius * 2,
        )
        span = (arc.end_angle - arc.start_angle) % 360
        if span == 0 and arc.end > arc.start:
            span = 359.99
        painter.drawArc(rect, int(-arc.start_angle * 16), int(-span * 16))

    def _draw_ticks(self, painter: object) -> None:
        painter.setPen(QPen(QColor("#212529"), 2))
        for tick in self._model.enzyme_ticks:
            painter.drawLine(_qpoint(tick.start_point), _qpoint(tick.end_point))

    def _draw_arrows(self, painter: object) -> None:
        for arrow in (*self._model.feature_arrows, *self._model.primer_arrows):
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QBrush(QColor(arrow.color)))
            painter.drawPolygon(QPolygonF([_qpoint(point) for point in arrow.points]))

    def _draw_scale_and_origin(self, painter: object) -> None:
        painter.setPen(QPen(QColor("#495057"), 2))
        painter.drawLine(_qpoint(self._model.scale.start_point), _qpoint(self._model.scale.end_point))
        painter.setPen(QPen(QColor("#d9480f"), 3))
        origin = self._model.origin_marker
        painter.drawLine(_qpoint(origin.start_point), _qpoint(origin.end_point))

    def _draw_labels(self, painter: object) -> None:
        painter.setPen(QPen(QColor("#212529"), 1))
        painter.setFont(QFont("Arial", 10))
        for label in self._model.labels:
            painter.drawText(QPointF(label.position.x, label.position.y), label.text)


def sample_plasmid_record() -> SequenceRecord:
    """Return a synthetic circular sample plasmid for the widget demo."""

    sequence = "ATGCGTAC" * 40
    return SequenceRecord(
        id="pSample",
        name="pSample",
        sequence=sequence,
        topology=MoleculeTopology.CIRCULAR,
        features=(
            Feature(type="promoter", start=12, end=48, strand=1, name="P_syn"),
            Feature(type="CDS", start=60, end=168, strand=1, name="synCds"),
            Feature(type="terminator", start=180, end=220, strand=1, name="T_syn"),
            Feature(type="rep_origin", start=245, end=310, strand=-1, name="ori"),
        ),
        primers=(
            Primer(name="Fwd", sequence="ATGCGTAC", start=60, end=68, strand=1),
            Primer(name="Rev", sequence="GTACGCAT", start=152, end=160, strand=-1),
        ),
        enzyme_sites=(
            EnzymeSite("EcoRI", "GAATTC", 32, 38, cut_index=1),
            EnzymeSite("BamHI", "GGATCC", 210, 216, cut_index=1),
        ),
    )


def show_sample_plasmid() -> int:
    """Open a minimal window showing a synthetic circular plasmid map."""

    if QApplication is None:
        msg = "PySide6 is required to open the sample plasmid widget"
        raise ImportError(msg)
    app = QApplication.instance() or QApplication(sys.argv)
    widget = PlasmidMapWidget(sample_plasmid_record())
    widget.show()
    return app.exec()


def _qpoint(point: Point) -> object:
    return QPointF(point.x, point.y)


if __name__ == "__main__":
    raise SystemExit(show_sample_plasmid())
