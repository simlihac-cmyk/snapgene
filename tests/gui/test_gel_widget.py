import os

import pytest

from plasmidlab.core import GelLaneInput, simulate_gel


def test_agarose_gel_widget_exports_svg_and_png(tmp_path) -> None:
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PySide6")

    from PySide6.QtWidgets import QApplication

    from plasmidlab.gui.gel_widget import AgaroseGelWidget

    model = simulate_gel(
        (GelLaneInput(name="digest", fragments=(3000, 1000, 500)),),
        ladder="1kb_dna_ladder",
    )
    svg_path = tmp_path / "gel.svg"
    png_path = tmp_path / "gel.png"

    app = QApplication.instance() or QApplication([])
    widget = AgaroseGelWidget(model)
    widget.write_svg(svg_path)
    widget.write_png(png_path)

    assert app is not None
    assert widget.model.lanes[0].is_ladder
    assert svg_path.read_text(encoding="utf-8").startswith('<?xml version="1.0"')
    assert png_path.read_bytes().startswith(b"\x89PNG")
