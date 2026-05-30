import os

import pytest

from plasmidlab.core import Feature, SequenceRecord


def test_gui_undo_redo_restores_exact_sequence_and_annotations() -> None:
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PySide6")

    from PySide6.QtWidgets import QApplication

    from plasmidlab.gui.main_window import PlasmidLabMainWindow

    original = SequenceRecord(
        id="undo",
        sequence="AAAA",
        features=(Feature(type="misc_feature", start=1, end=3, name="original"),),
    )
    edited = original.insert(2, "GG")

    app = QApplication.instance() or QApplication([])
    window = PlasmidLabMainWindow((original,))
    window._set_current_record(edited)

    assert app is not None
    assert window.current_record == edited
    assert window.undo_action.isEnabled()
    assert not window.redo_action.isEnabled()

    window.undo()

    assert window.current_record == original
    assert window.current_record.features == original.features
    assert not window.undo_action.isEnabled()
    assert window.redo_action.isEnabled()

    window.redo()

    assert window.current_record == edited
    assert window.current_record.features == edited.features
    assert window.undo_action.isEnabled()
    assert not window.redo_action.isEnabled()
