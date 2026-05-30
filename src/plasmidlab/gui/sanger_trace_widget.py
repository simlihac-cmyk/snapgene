"""PySide6 chromatogram widget for parsed Sanger traces."""

from __future__ import annotations

import sys
from pathlib import Path

from plasmidlab.io.sanger import SangerTrace, read_ab1

try:
    from PySide6.QtCore import QPointF, Qt
    from PySide6.QtGui import QColor, QFont, QPainter, QPen, QPolygonF
    from PySide6.QtWidgets import QApplication, QMainWindow, QWidget
except ImportError:  # pragma: no cover - exercised only without PySide6.
    QApplication = None  # type: ignore[assignment]
    QWidget = object  # type: ignore[assignment,misc]
    QMainWindow = object  # type: ignore[assignment,misc]
    QPointF = None  # type: ignore[assignment]
    Qt = None  # type: ignore[assignment]
    QColor = None  # type: ignore[assignment]
    QFont = None  # type: ignore[assignment]
    QPainter = None  # type: ignore[assignment]
    QPen = None  # type: ignore[assignment]
    QPolygonF = None  # type: ignore[assignment]


TRACE_COLORS = {
    "A": "#2f9e44",
    "C": "#1971c2",
    "G": "#212529",
    "T": "#d9480f",
}


class SangerTraceWidget(QWidget):
    """Render called bases and chromatogram channels from a Sanger trace."""

    def __init__(self, trace: SangerTrace | None = None, parent: object | None = None) -> None:
        if QApplication is None:
            msg = "PySide6 is required to use SangerTraceWidget"
            raise ImportError(msg)
        super().__init__(parent)  # type: ignore[misc]
        self._trace = trace or sample_trace()
        self.setMinimumSize(900, 320)
        self.setWindowTitle(f"{self._trace.id} Sanger trace")

    @property
    def trace(self) -> SangerTrace:
        """Return the trace currently rendered by this widget."""

        return self._trace

    def set_trace(self, trace: SangerTrace) -> None:
        """Replace the displayed trace and schedule repaint."""

        self._trace = trace
        self.setWindowTitle(f"{trace.id} Sanger trace")
        self.update()

    def paintEvent(self, event: object) -> None:  # noqa: N802 - Qt API name.
        del event
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.fillRect(self.rect(), QColor("#ffffff"))
        self._draw_axes(painter)
        self._draw_channels(painter)
        self._draw_base_calls(painter)
        painter.end()

    def _draw_axes(self, painter: object) -> None:
        rect = self.rect()
        left = 42
        right = rect.width() - 18
        baseline = rect.height() - 92
        painter.setPen(QPen(QColor("#adb5bd"), 1))
        painter.drawLine(QPointF(left, baseline), QPointF(right, baseline))

    def _draw_channels(self, painter: object) -> None:
        if not self._trace.chromatogram:
            return
        rect = self.rect()
        left = 42.0
        right = float(rect.width() - 18)
        top = 28.0
        baseline = float(rect.height() - 92)
        channel_length = max((len(values) for values in self._trace.chromatogram.values()), default=0)
        max_signal = max((max(values) for values in self._trace.chromatogram.values() if values), default=1)
        if channel_length < 2 or max_signal <= 0:
            return
        x_scale = (right - left) / (channel_length - 1)
        y_scale = (baseline - top) / max_signal
        for base in ("G", "A", "T", "C"):
            values = self._trace.chromatogram.get(base)
            if not values:
                continue
            painter.setPen(QPen(QColor(TRACE_COLORS[base]), 1.4))
            points = [
                QPointF(left + index * x_scale, baseline - signal * y_scale)
                for index, signal in enumerate(values)
            ]
            painter.drawPolyline(QPolygonF(points))

    def _draw_base_calls(self, painter: object) -> None:
        rect = self.rect()
        left = 42.0
        right = float(rect.width() - 18)
        base_y = float(rect.height() - 64)
        quality_y = float(rect.height() - 34)
        channel_length = max((len(values) for values in self._trace.chromatogram.values()), default=0)
        if channel_length > 1:
            def x_for_peak(peak: int) -> float:
                return left + (right - left) * peak / (channel_length - 1)
        else:
            spacing = (right - left) / max(len(self._trace.called_bases), 1)

            def x_for_peak(peak: int) -> float:
                return left + spacing * peak

        painter.setFont(QFont("Courier New", 11))
        for index, base in enumerate(self._trace.called_bases):
            peak = self._trace.peak_positions[index] if index < len(self._trace.peak_positions) else index
            x = x_for_peak(peak)
            painter.setPen(QPen(QColor(TRACE_COLORS.get(base, "#495057")), 1))
            painter.drawText(QPointF(x - 4, base_y), base)
            if index < len(self._trace.qualities):
                painter.setPen(QPen(QColor("#868e96"), 1))
                painter.drawText(QPointF(x - 8, quality_y), str(self._trace.qualities[index]))


class SangerTraceWindow(QMainWindow):
    """Small standalone window for a Sanger trace widget."""

    def __init__(self, trace: SangerTrace, parent: object | None = None) -> None:
        if QApplication is None:
            msg = "PySide6 is required to use SangerTraceWindow"
            raise ImportError(msg)
        super().__init__(parent)  # type: ignore[misc]
        self.trace_widget = SangerTraceWidget(trace)
        self.setCentralWidget(self.trace_widget)
        self.resize(980, 380)
        self.setWindowTitle(f"{trace.id} Sanger trace")


def sample_trace() -> SangerTrace:
    """Return a deterministic synthetic trace for smoke tests and demos."""

    called_bases = "ACGTACGT"
    peak_positions = tuple(index * 18 + 8 for index in range(len(called_bases)))
    length = peak_positions[-1] + 22
    chromatogram: dict[str, tuple[int, ...]] = {}
    for base in "ACGT":
        values: list[int] = []
        for index in range(length):
            signal = 0
            for peak, called_base in zip(peak_positions, called_bases, strict=True):
                if called_base == base:
                    signal = max(signal, max(0, 120 - abs(index - peak) * 18))
            values.append(signal)
        chromatogram[base] = tuple(values)
    return SangerTrace(
        id="synthetic_trace",
        called_bases=called_bases,
        qualities=(35, 36, 34, 33, 32, 31, 30, 29),
        peak_positions=peak_positions,
        chromatogram=chromatogram,
    )


def open_trace_window(path: str | Path) -> int:
    """Open a Sanger trace file in a standalone chromatogram window."""

    if QApplication is None:
        msg = "PySide6 is required to open Sanger trace files"
        raise ImportError(msg)
    app = QApplication.instance() or QApplication(sys.argv)
    window = SangerTraceWindow(read_ab1(path))
    window.show()
    return app.exec()


if __name__ == "__main__":
    if len(sys.argv) > 1:
        raise SystemExit(open_trace_window(sys.argv[1]))
    app = QApplication.instance() or QApplication(sys.argv)
    widget = SangerTraceWidget(sample_trace())
    widget.show()
    raise SystemExit(app.exec())
