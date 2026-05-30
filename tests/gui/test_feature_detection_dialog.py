import os

import pytest

from plasmidlab.core import FeatureLibraryEntry, SequenceRecord, detect_features


def test_feature_detection_dialog_allows_review_before_apply() -> None:
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PySide6")

    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QApplication

    from plasmidlab.gui.main_window import FeatureDetectionDialog

    entry = FeatureLibraryEntry(name="review feature", type="misc_feature", sequence="ATGAAACCC")
    detections = detect_features(SequenceRecord(id="review", sequence="ATGAAACCC"), (entry,))

    app = QApplication.instance() or QApplication([])
    dialog = FeatureDetectionDialog(None, detections)
    dialog.table.item(0, 0).setCheckState(Qt.CheckState.Unchecked)

    assert app is not None
    assert dialog.table.rowCount() == 1
    assert dialog.selected_detections() == ()


def test_main_window_has_detect_features_menu_and_can_apply_proposals() -> None:
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PySide6")

    from PySide6.QtWidgets import QApplication

    from plasmidlab.gui.main_window import PlasmidLabMainWindow

    entry = FeatureLibraryEntry(name="gui feature", type="misc_feature", sequence="ATGAAACCC")
    record = SequenceRecord(id="gui_detect", sequence="GGGATGAAACCCTTT")
    detections = detect_features(record, (entry,))

    app = QApplication.instance() or QApplication([])
    window = PlasmidLabMainWindow((record,))
    window._apply_feature_detections(detections)

    assert app is not None
    assert window.detect_features_action.text() == "Detect Features"
    assert window.current_record is not None
    assert len(window.current_record.features) == 1
    assert window.feature_table.rowCount() == 1
    assert window.history_table.rowCount() == 1
