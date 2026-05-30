"""Simple PySide6 visualization for PlasmidLab history graphs."""

from __future__ import annotations

from plasmidlab.history import HistoryGraph, build_history_graph

try:
    from PySide6.QtCore import QPointF, QRectF, Qt
    from PySide6.QtGui import QColor, QFont, QPainter, QPen
    from PySide6.QtWidgets import QApplication, QWidget
except ImportError:  # pragma: no cover - exercised only without PySide6.
    QApplication = None  # type: ignore[assignment]
    QWidget = object  # type: ignore[assignment,misc]
    QPointF = None  # type: ignore[assignment]
    QRectF = None  # type: ignore[assignment]
    Qt = None  # type: ignore[assignment]
    QColor = None  # type: ignore[assignment]
    QFont = None  # type: ignore[assignment]
    QPainter = None  # type: ignore[assignment]
    QPen = None  # type: ignore[assignment]


class HistoryGraphWidget(QWidget):
    """Draw record and operation nodes with derivation edges."""

    def __init__(self, graph: HistoryGraph | None = None, parent: object | None = None) -> None:
        if QApplication is None:
            msg = "PySide6 is required to use HistoryGraphWidget"
            raise ImportError(msg)
        super().__init__(parent)  # type: ignore[misc]
        self._graph = graph or build_history_graph(())
        self.setMinimumHeight(180)

    @property
    def graph(self) -> HistoryGraph:
        """Return the displayed graph."""

        return self._graph

    def set_graph(self, graph: HistoryGraph) -> None:
        """Replace the graph and schedule repaint."""

        self._graph = graph
        self.update()

    def paintEvent(self, event: object) -> None:  # noqa: N802 - Qt API name.
        del event
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.fillRect(self.rect(), QColor("#ffffff"))
        positions = self._node_positions()
        painter.setFont(QFont("Arial", 9))
        painter.setPen(QPen(QColor("#adb5bd"), 1.5))
        for edge in self._graph.edges:
            if edge.source in positions and edge.target in positions:
                painter.drawLine(_center(positions[edge.source]), _center(positions[edge.target]))
        for node in self._graph.nodes:
            rect = positions.get(node.id)
            if rect is None:
                continue
            fill = QColor("#e7f5ff" if node.kind == "record" else "#fff3bf")
            stroke = QColor("#1c7ed6" if node.kind == "record" else "#f08c00")
            painter.setPen(QPen(stroke, 1.5))
            painter.setBrush(fill)
            painter.drawRoundedRect(rect, 5, 5)
            painter.setPen(QPen(QColor("#212529"), 1))
            painter.drawText(rect.adjusted(5, 2, -5, -2), Qt.AlignmentFlag.AlignCenter, node.label)
        painter.end()

    def _node_positions(self) -> dict[str, object]:
        width = max(320, self.width())
        y_record = 26
        y_event = 102
        record_nodes = [node for node in self._graph.nodes if node.kind == "record"]
        event_nodes = [node for node in self._graph.nodes if node.kind == "event"]
        positions: dict[str, object] = {}
        for index, node in enumerate(record_nodes):
            x = _row_x(index, len(record_nodes), width)
            positions[node.id] = QRectF(x - 60, y_record, 120, 34)
        for index, node in enumerate(event_nodes):
            x = _row_x(index, len(event_nodes), width)
            positions[node.id] = QRectF(x - 72, y_event, 144, 40)
        return positions


def _row_x(index: int, count: int, width: int) -> float:
    if count <= 1:
        return width / 2
    margin = 90
    return margin + (width - margin * 2) * index / (count - 1)


def _center(rect: object) -> object:
    return QPointF(rect.x() + rect.width() / 2, rect.y() + rect.height() / 2)
