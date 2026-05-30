import os

import pytest

from plasmidlab.io.sanger import SangerTrace


def test_sanger_trace_widget_renders_mock_trace() -> None:
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PySide6")

    from PySide6.QtGui import QImage
    from PySide6.QtWidgets import QApplication

    from plasmidlab.gui.sanger_trace_widget import SangerTraceWidget

    trace = SangerTrace(
        id="mock",
        called_bases="ACGT",
        qualities=(30, 31, 32, 33),
        peak_positions=(4, 12, 20, 28),
        chromatogram={
            "A": (0, 20, 80, 20, 0),
            "C": (0, 0, 20, 80, 20),
            "G": (20, 80, 20, 0, 0),
            "T": (0, 0, 0, 20, 80),
        },
    )

    app = QApplication.instance() or QApplication([])
    widget = SangerTraceWidget(trace)
    image = QImage(900, 320, QImage.Format.Format_ARGB32)
    widget.resize(900, 320)
    widget.render(image)

    assert app is not None
    assert widget.trace.called_bases == "ACGT"
    assert widget.trace.qualities == (30, 31, 32, 33)
