import os

import pytest

from plasmidlab.core import Feature, FeatureSegment, MoleculeTopology, SequenceRecord
from plasmidlab.io.genbank import write_genbank


def test_main_window_opens_genbank_shows_tabs_and_runs_enzyme_search(tmp_path) -> None:
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PySide6")

    from PySide6.QtWidgets import QApplication

    from plasmidlab.gui.main_window import PlasmidLabMainWindow

    sequence = "A" * 10 + "GAATTC" + "C" * 20 + "GGATCC" + "T" * 20
    record = SequenceRecord(
        id="gui_record",
        sequence=sequence,
        topology=MoleculeTopology.CIRCULAR,
        features=(Feature(type="CDS", start=0, end=9, strand=1, name="tiny_cds"),),
    )
    path = tmp_path / "record.gb"
    svg_path = tmp_path / "record.svg"
    gel_svg_path = tmp_path / "gel.svg"
    gel_png_path = tmp_path / "gel.png"
    history_path = tmp_path / "history.json"
    saved_path = tmp_path / "saved.gb"
    write_genbank(record, path)

    app = QApplication.instance() or QApplication([])
    window = PlasmidLabMainWindow()
    window.open_path(path)
    window.project_search_edit.setText("tiny_cds")
    window.search_project_panel()
    window.enzyme_combo.setCurrentText("EcoRI,BamHI,HindIII")
    window.run_enzyme_search()
    window.export_svg_path(svg_path)
    window.export_gel_svg_path(gel_svg_path)
    window.export_gel_png_path(gel_png_path)
    window.export_history_json_path(history_path)
    window.save_path(saved_path)

    assert app is not None
    assert window.current_record is not None
    assert window.current_record.id == "gui_record"
    assert window.tabs.count() == 7
    assert window.project_list.count() == 1
    assert window.project_search_edit.placeholderText() == "Search project"
    assert window.map_widget.model.map_kind == "circular"
    assert "GAATTC" in window.sequence_view.toPlainText()
    assert window.feature_table.rowCount() == 1
    assert window.enzyme_table.rowCount() == 2
    assert "HindIII" in window.enzyme_summary.text()
    assert window.gel_widget.model.lanes[0].is_ladder
    assert window.gel_widget.model.lanes[1].name == "Digest"
    assert svg_path.read_text(encoding="utf-8").startswith('<?xml version="1.0"')
    assert gel_svg_path.read_text(encoding="utf-8").startswith('<?xml version="1.0"')
    assert gel_png_path.read_bytes().startswith(b"\x89PNG")
    assert '"nodes"' in history_path.read_text(encoding="utf-8")
    assert saved_path.exists()


def test_main_window_does_not_simplify_compound_feature_in_simple_editor(monkeypatch) -> None:
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PySide6")

    from PySide6.QtWidgets import QApplication, QMessageBox

    from plasmidlab.gui.main_window import PlasmidLabMainWindow

    record = SequenceRecord(
        id="compound_gui",
        sequence="A" * 100,
        topology=MoleculeTopology.CIRCULAR,
        features=(
            Feature(
                type="CDS",
                name="split",
                segments=(FeatureSegment(80, 100, strand=1), FeatureSegment(0, 10, strand=1)),
            ),
        ),
    )
    messages: list[str] = []

    def fake_information(parent: object, title: str, text: str) -> object:
        _ = (parent, title)
        messages.append(text)
        return QMessageBox.StandardButton.Ok

    monkeypatch.setattr(QMessageBox, "information", fake_information)
    app = QApplication.instance() or QApplication([])
    window = PlasmidLabMainWindow()
    window.load_records((record,))
    window.feature_table.setCurrentCell(0, 0)

    window.edit_feature()

    assert app is not None
    assert window.current_record is not None
    assert window.current_record.features[0].segments == record.features[0].segments
    assert messages == ["Compound features require advanced editing and are not editable in this dialog yet."]
