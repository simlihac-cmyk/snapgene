import os

import pytest


def test_plasmid_map_widget_constructs_sample_record() -> None:
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PySide6")

    from PySide6.QtWidgets import QApplication

    from plasmidlab.gui.plasmid_map_widget import PlasmidMapWidget, sample_plasmid_record

    app = QApplication.instance() or QApplication([])
    widget = PlasmidMapWidget(sample_plasmid_record())

    assert app is not None
    assert widget.model.map_kind == "circular"
    assert widget.model.record_id == "pSample"
    assert widget.model.feature_arcs
