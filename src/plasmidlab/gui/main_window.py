"""First usable PySide6 desktop shell for PlasmidLab."""

from __future__ import annotations

import json
import sys
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from plasmidlab.core import (
    Feature,
    FeatureDetection,
    GelInputFragment,
    GelLaneInput,
    Primer,
    ProvenanceEvent,
    SequenceRecord,
    analyze_restriction_sites,
    apply_feature_annotations,
    detect_features,
    digest,
    find_primer_bindings,
    primer_metrics,
    simulate_gel,
    translate_feature,
)
from plasmidlab.core.restriction import RestrictionAnalysis
from plasmidlab.gui.gel_widget import AgaroseGelWidget
from plasmidlab.gui.history_graph_widget import HistoryGraphWidget
from plasmidlab.gui.plasmid_map_widget import PlasmidMapWidget
from plasmidlab.gui.sanger_trace_widget import SangerTraceWindow
from plasmidlab.history import build_history_graph, write_history_json
from plasmidlab.io.fasta import read_fasta, write_fasta
from plasmidlab.io.genbank import read_genbank, write_genbank
from plasmidlab.io.sanger import read_ab1
from plasmidlab.project import ProjectDatabase
from plasmidlab.render import prepare_map_overlays, render_plasmid_map, write_svg

try:
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QAction, QColor, QFont, QTextCharFormat, QTextCursor
    from PySide6.QtWidgets import (
        QApplication,
        QComboBox,
        QDialog,
        QDialogButtonBox,
        QFileDialog,
        QFormLayout,
        QHBoxLayout,
        QHeaderView,
        QLabel,
        QLineEdit,
        QListWidget,
        QMainWindow,
        QMessageBox,
        QPushButton,
        QPlainTextEdit,
        QSpinBox,
        QSplitter,
        QTableWidget,
        QTableWidgetItem,
        QTabWidget,
        QTextEdit,
        QVBoxLayout,
        QWidget,
    )
except ImportError:  # pragma: no cover - lets the package import without GUI deps.
    QApplication = None  # type: ignore[assignment]
    QMainWindow = object  # type: ignore[assignment,misc]
    QDialog = object  # type: ignore[assignment,misc]
    QWidget = object  # type: ignore[assignment,misc]

    def QColor(value: str) -> str:  # type: ignore[no-redef]
        return value


COMMON_ENZYME_SETS = {
    "Common": "EcoRI,BamHI,HindIII,XhoI,NotI,SalI,KpnI,SmaI,PstI",
    "EcoRI/BamHI": "EcoRI,BamHI",
    "Cloning trio": "EcoRI,BamHI,HindIII",
}

FEATURE_COLORS = {
    "promoter": QColor("#d3f9d8"),
    "CDS": QColor("#d0ebff"),
    "terminator": QColor("#ffe3e3"),
    "rep_origin": QColor("#e5dbff"),
    "misc_feature": QColor("#f1f3f5"),
}


class PlasmidLabMainWindow(QMainWindow):
    """Main desktop window backed by PlasmidLab core APIs."""

    def __init__(self, records: tuple[SequenceRecord, ...] | None = None) -> None:
        if QApplication is None:
            msg = "PySide6 is required to use PlasmidLabMainWindow"
            raise ImportError(msg)
        super().__init__()
        self._records: list[SequenceRecord] = []
        self._current_index = -1
        self._current_path: Path | None = None
        self._last_restriction_analysis: RestrictionAnalysis | None = None
        self._last_primer_bindings: dict[int, tuple[Any, ...]] = {}
        self._last_digest_fragments: tuple[Any, ...] | None = None
        self._trace_windows: list[SangerTraceWindow] = []
        self._undo_stack: list[SequenceRecord] = []
        self._redo_stack: list[SequenceRecord] = []
        self._project_db: ProjectDatabase | None = None
        self._project_record_pks: list[int] = []
        self._project_list_indices: list[int] = []

        self.setWindowTitle("PlasmidLab")
        self.resize(1280, 820)
        self._build_actions()
        self._build_ui()
        if records:
            self.load_records(records)

    @property
    def current_record(self) -> SequenceRecord | None:
        """Return the selected sequence record."""

        if self._current_index < 0 or self._current_index >= len(self._records):
            return None
        return self._records[self._current_index]

    def load_records(
        self,
        records: tuple[SequenceRecord, ...],
        *,
        source: Path | None = None,
        project_pks: tuple[int, ...] | None = None,
    ) -> None:
        """Load records into the project panel and show the first record."""

        self._records = list(records)
        self._project_record_pks = list(project_pks or ())
        self._current_path = source
        self._current_index = 0 if self._records else -1
        self._last_restriction_analysis = None
        self._last_digest_fragments = None
        self._last_primer_bindings = {}
        self._undo_stack.clear()
        self._redo_stack.clear()
        self._populate_project_list()
        if self._records:
            self.project_list.setCurrentRow(0)
        self._refresh_all()
        self._update_undo_redo_actions()

    def open_path(self, path: str | Path) -> None:
        """Open a FASTA, GenBank, or AB1 file without showing a file dialog."""

        input_path = Path(path)
        if input_path.suffix.lower() in {".ab1", ".abi"}:
            self.open_trace_path(input_path)
            return
        if input_path.suffix.lower() in {".fa", ".fasta", ".fna"}:
            records = read_fasta(input_path)
        else:
            records = read_genbank(input_path)
        self.load_records(records, source=input_path)

    def open_trace_path(self, path: str | Path) -> SangerTraceWindow:
        """Open an AB1 trace in a chromatogram window without base calling."""

        trace_window = SangerTraceWindow(read_ab1(path), parent=self)
        self._trace_windows.append(trace_window)
        trace_window.show()
        return trace_window

    def open_project_database_path(self, path: str | Path) -> None:
        """Open a SQLite project database without showing a file dialog."""

        if self._project_db is not None:
            self._project_db.close()
        self._project_db = ProjectDatabase(path)
        self.load_records(
            self._project_db.all_records(),
            source=Path(path),
            project_pks=self._project_db.record_pks(),
        )

    def import_folder_to_project_path(self, folder: str | Path) -> tuple[int, ...]:
        """Import FASTA/GenBank files into the open project database."""

        if self._project_db is None:
            self._project_db = ProjectDatabase(":memory:")
        imported = self._project_db.import_folder(folder, replace=True)
        self.load_records(
            self._project_db.all_records(),
            source=Path(folder),
            project_pks=self._project_db.record_pks(),
        )
        return imported

    def save_path(self, path: str | Path) -> None:
        """Save the current record as FASTA or GenBank without showing a file dialog."""

        record = self.current_record
        if record is None:
            return
        output_path = Path(path)
        if output_path.suffix.lower() in {".fa", ".fasta", ".fna"}:
            write_fasta(record, output_path)
        else:
            write_genbank(record, output_path)

    def export_svg_path(self, path: str | Path) -> None:
        """Export the current map SVG without showing a file dialog."""

        record = self.current_record
        if record is None:
            return
        write_svg(
            render_plasmid_map(
                record,
                overlays=self._map_overlays(record),
            ),
            path,
        )

    def export_gel_svg_path(self, path: str | Path) -> None:
        """Export the current gel SVG without showing a file dialog."""

        self.gel_widget.write_svg(path)

    def export_gel_png_path(self, path: str | Path) -> None:
        """Export the current gel PNG without showing a file dialog."""

        self.gel_widget.write_png(path)

    def export_history_json_path(self, path: str | Path) -> None:
        """Export project history as JSON without showing a file dialog."""

        write_history_json(self._records, path)

    def run_enzyme_search(self) -> None:
        """Run restriction enzyme search using the selected enzyme text."""

        record = self.current_record
        if record is None:
            return
        enzyme_text = self.enzyme_combo.currentText().strip()
        if not enzyme_text:
            return
        try:
            analysis = analyze_restriction_sites(record, enzyme_text)
        except Exception as error:  # pragma: no cover - dialog path.
            self._show_error("Enzyme Search Failed", str(error))
            return
        self._last_restriction_analysis = analysis
        self._last_digest_fragments = digest(record, enzyme_text)
        self._populate_enzyme_table(analysis)
        self._refresh_map()
        self._refresh_gel()

    def _build_actions(self) -> None:
        file_menu = self.menuBar().addMenu("&File")
        open_action = QAction("&Open FASTA/GenBank/AB1...", self)
        open_action.triggered.connect(self.open_file)
        save_action = QAction("&Save As...", self)
        save_action.triggered.connect(self.save_file)
        export_svg_action = QAction("Export Map &SVG...", self)
        export_svg_action.triggered.connect(self.export_svg)
        export_gel_svg_action = QAction("Export Gel SVG...", self)
        export_gel_svg_action.triggered.connect(self.export_gel_svg)
        export_gel_png_action = QAction("Export Gel PNG...", self)
        export_gel_png_action.triggered.connect(self.export_gel_png)
        export_history_action = QAction("Export History JSON...", self)
        export_history_action.triggered.connect(self.export_history_json)
        exit_action = QAction("E&xit", self)
        exit_action.triggered.connect(self.close)
        file_menu.addAction(open_action)
        open_project_action = QAction("Open Project Database...", self)
        open_project_action.triggered.connect(self.open_project_database)
        import_folder_action = QAction("Import Folder to Project...", self)
        import_folder_action.triggered.connect(self.import_folder_to_project)
        file_menu.addAction(open_project_action)
        file_menu.addAction(import_folder_action)
        file_menu.addAction(save_action)
        file_menu.addAction(export_svg_action)
        file_menu.addAction(export_gel_svg_action)
        file_menu.addAction(export_gel_png_action)
        file_menu.addAction(export_history_action)
        file_menu.addSeparator()
        file_menu.addAction(exit_action)

        edit_menu = self.menuBar().addMenu("&Edit")
        self.undo_action = QAction("&Undo", self)
        self.undo_action.triggered.connect(self.undo)
        self.redo_action = QAction("&Redo", self)
        self.redo_action.triggered.connect(self.redo)
        edit_menu.addAction(self.undo_action)
        edit_menu.addAction(self.redo_action)
        self._update_undo_redo_actions()

        tools_menu = self.menuBar().addMenu("&Tools")
        enzyme_action = QAction("Run &Enzyme Search", self)
        enzyme_action.triggered.connect(self.run_enzyme_search)
        self.detect_features_action = QAction("Detect Features", self)
        self.detect_features_action.triggered.connect(self.run_feature_detection)
        tools_menu.addAction(enzyme_action)
        tools_menu.addAction(self.detect_features_action)

    def _build_ui(self) -> None:
        splitter = QSplitter(Qt.Orientation.Horizontal)
        project_panel = QWidget()
        project_layout = QVBoxLayout(project_panel)
        search_row = QHBoxLayout()
        self.project_search_edit = QLineEdit()
        self.project_search_edit.setPlaceholderText("Search project")
        self.project_search_edit.returnPressed.connect(self.search_project_panel)
        search_button = QPushButton("Search")
        search_button.clicked.connect(self.search_project_panel)
        clear_button = QPushButton("Clear")
        clear_button.clicked.connect(self.clear_project_search)
        search_row.addWidget(self.project_search_edit, 1)
        search_row.addWidget(search_button)
        search_row.addWidget(clear_button)
        self.project_list = QListWidget()
        self.project_list.setMinimumWidth(230)
        self.project_list.currentRowChanged.connect(self._record_selection_changed)
        self.project_status_label = QLabel("No project database open.")
        project_layout.addLayout(search_row)
        project_layout.addWidget(self.project_list, 1)
        project_layout.addWidget(self.project_status_label)
        splitter.addWidget(project_panel)

        self.tabs = QTabWidget()
        self.map_widget = PlasmidMapWidget()
        self.tabs.addTab(self.map_widget, "Map")
        self.tabs.addTab(self._build_sequence_tab(), "Sequence")
        self.tabs.addTab(self._build_features_tab(), "Features")
        self.tabs.addTab(self._build_enzymes_tab(), "Enzymes")
        self.tabs.addTab(self._build_primers_tab(), "Primers")
        self.gel_widget = AgaroseGelWidget()
        self.tabs.addTab(self.gel_widget, "Gel")
        self.tabs.addTab(self._build_history_tab(), "History")
        splitter.addWidget(self.tabs)
        splitter.setStretchFactor(1, 1)
        self.setCentralWidget(splitter)
        self.statusBar().showMessage("Ready")

    def _build_sequence_tab(self) -> QWidget:
        widget = QWidget()
        layout = QVBoxLayout(widget)
        self.sequence_view = QPlainTextEdit()
        self.sequence_view.setReadOnly(True)
        self.sequence_view.setFont(QFont("Courier New", 10))
        self.sequence_view.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.sequence_view.selectionChanged.connect(self._update_selection_status)
        self.translation_view = QPlainTextEdit()
        self.translation_view.setReadOnly(True)
        self.translation_view.setFont(QFont("Courier New", 10))
        self.translation_view.setMaximumHeight(150)
        self.selection_label = QLabel("Selection: none")
        layout.addWidget(self.selection_label)
        layout.addWidget(self.sequence_view, 3)
        layout.addWidget(QLabel("CDS translations"))
        layout.addWidget(self.translation_view, 1)
        return widget

    def _build_features_tab(self) -> QWidget:
        widget = QWidget()
        layout = QVBoxLayout(widget)
        button_row = QHBoxLayout()
        add_button = QPushButton("Add")
        edit_button = QPushButton("Edit")
        delete_button = QPushButton("Delete")
        add_button.clicked.connect(self.add_feature)
        edit_button.clicked.connect(self.edit_feature)
        delete_button.clicked.connect(self.delete_selected_features)
        button_row.addWidget(add_button)
        button_row.addWidget(edit_button)
        button_row.addWidget(delete_button)
        button_row.addStretch()
        self.feature_table = QTableWidget(0, 6)
        self.feature_table.setHorizontalHeaderLabels(["Name", "Type", "Strand", "Start", "End", "Qualifiers"])
        self.feature_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.feature_table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        layout.addLayout(button_row)
        layout.addWidget(self.feature_table)
        return widget

    def _build_enzymes_tab(self) -> QWidget:
        widget = QWidget()
        layout = QVBoxLayout(widget)
        control_row = QHBoxLayout()
        self.enzyme_combo = QComboBox()
        self.enzyme_combo.setEditable(True)
        for label, enzymes in COMMON_ENZYME_SETS.items():
            self.enzyme_combo.addItem(enzymes, label)
        search_button = QPushButton("Analyze")
        search_button.clicked.connect(self.run_enzyme_search)
        control_row.addWidget(QLabel("Enzymes"))
        control_row.addWidget(self.enzyme_combo, 1)
        control_row.addWidget(search_button)
        self.enzyme_summary = QLabel("No enzyme search has been run.")
        self.enzyme_table = QTableWidget(0, 7)
        self.enzyme_table.setHorizontalHeaderLabels(
            ["Enzyme", "Start", "End", "Strand", "Top cut", "Bottom cut", "Recognition"]
        )
        self.enzyme_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addLayout(control_row)
        layout.addWidget(self.enzyme_summary)
        layout.addWidget(self.enzyme_table)
        return widget

    def _build_primers_tab(self) -> QWidget:
        widget = QWidget()
        layout = QVBoxLayout(widget)
        button_row = QHBoxLayout()
        add_button = QPushButton("Add Primer")
        refresh_button = QPushButton("Refresh Bindings")
        add_button.clicked.connect(self.add_primer)
        refresh_button.clicked.connect(self._refresh_primers)
        button_row.addWidget(add_button)
        button_row.addWidget(refresh_button)
        button_row.addStretch()
        self.primer_table = QTableWidget(0, 7)
        self.primer_table.setHorizontalHeaderLabels(
            ["Name", "Sequence", "Length", "GC %", "Wallace Tm", "Bindings", "Binding Sites"]
        )
        self.primer_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addLayout(button_row)
        layout.addWidget(self.primer_table)
        return widget

    def _build_history_tab(self) -> QWidget:
        widget = QWidget()
        layout = QVBoxLayout(widget)
        splitter = QSplitter(Qt.Orientation.Vertical)
        self.history_graph_widget = HistoryGraphWidget()
        self.history_table = QTableWidget(0, 6)
        self.history_table.setHorizontalHeaderLabels(
            ["Time", "Operation", "Inputs", "Output", "Ranges", "Description"]
        )
        self.history_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        splitter.addWidget(self.history_graph_widget)
        splitter.addWidget(self.history_table)
        splitter.setStretchFactor(1, 1)
        layout.addWidget(splitter)
        return widget

    def open_file(self) -> None:
        path_text, _ = QFileDialog.getOpenFileName(
            self,
            "Open Sequence",
            "",
            "Sequence and trace files (*.gb *.gbk *.genbank *.fa *.fasta *.fna *.ab1 *.abi);;All files (*)",
        )
        if not path_text:
            return
        path = Path(path_text)
        try:
            self.open_path(path)
        except Exception as error:  # pragma: no cover - dialog path.
            self._show_error("Open Failed", str(error))

    def open_project_database(self) -> None:
        path_text, _ = QFileDialog.getOpenFileName(
            self,
            "Open Project Database",
            "",
            "SQLite databases (*.sqlite *.sqlite3 *.db);;All files (*)",
        )
        if not path_text:
            return
        try:
            self.open_project_database_path(path_text)
        except Exception as error:  # pragma: no cover - dialog path.
            self._show_error("Open Project Failed", str(error))

    def import_folder_to_project(self) -> None:
        folder_text = QFileDialog.getExistingDirectory(self, "Import Folder to Project")
        if not folder_text:
            return
        try:
            imported = self.import_folder_to_project_path(folder_text)
        except Exception as error:  # pragma: no cover - dialog path.
            self._show_error("Import Failed", str(error))
            return
        self.statusBar().showMessage(f"Imported {len(imported)} records")

    def save_file(self) -> None:
        record = self.current_record
        if record is None:
            return
        path_text, _ = QFileDialog.getSaveFileName(
            self,
            "Save Sequence",
            str(self._current_path or Path(f"{record.id}.gb")),
            "GenBank (*.gb);;FASTA (*.fasta)",
        )
        if not path_text:
            return
        try:
            self.save_path(path_text)
        except Exception as error:  # pragma: no cover - dialog path.
            self._show_error("Save Failed", str(error))

    def export_svg(self) -> None:
        record = self.current_record
        if record is None:
            return
        path_text, _ = QFileDialog.getSaveFileName(
            self,
            "Export SVG",
            f"{record.id}.svg",
            "SVG (*.svg)",
        )
        if not path_text:
            return
        try:
            self.export_svg_path(path_text)
        except Exception as error:  # pragma: no cover - dialog path.
            self._show_error("Export Failed", str(error))

    def export_gel_svg(self) -> None:
        path_text, _ = QFileDialog.getSaveFileName(
            self,
            "Export Gel SVG",
            "agarose_gel.svg",
            "SVG (*.svg)",
        )
        if not path_text:
            return
        try:
            self.export_gel_svg_path(path_text)
        except Exception as error:  # pragma: no cover - dialog path.
            self._show_error("Export Failed", str(error))

    def export_gel_png(self) -> None:
        path_text, _ = QFileDialog.getSaveFileName(
            self,
            "Export Gel PNG",
            "agarose_gel.png",
            "PNG (*.png)",
        )
        if not path_text:
            return
        try:
            self.export_gel_png_path(path_text)
        except Exception as error:  # pragma: no cover - dialog path.
            self._show_error("Export Failed", str(error))

    def export_history_json(self) -> None:
        path_text, _ = QFileDialog.getSaveFileName(
            self,
            "Export History JSON",
            "plasmidlab_history.json",
            "JSON (*.json)",
        )
        if not path_text:
            return
        try:
            self.export_history_json_path(path_text)
        except Exception as error:  # pragma: no cover - dialog path.
            self._show_error("Export Failed", str(error))

    def add_feature(self) -> None:
        record = self.current_record
        if record is None:
            return
        dialog = FeatureDialog(self, sequence_length=record.length)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            feature = dialog.feature()
            self._replace_current_record(
                features=record.features + (feature,),
                operation="add_feature",
                parameters={"name": feature.name or "", "type": feature.type},
            )

    def edit_feature(self) -> None:
        record = self.current_record
        if record is None:
            return
        row = self.feature_table.currentRow()
        if row < 0 or row >= len(record.features):
            return
        if not _feature_editable_in_simple_dialog(record.features[row]):
            QMessageBox.information(
                self,
                "Feature Editing",
                "Compound features require advanced editing and are not editable in this dialog yet.",
            )
            return
        dialog = FeatureDialog(self, sequence_length=record.length, feature=record.features[row])
        if dialog.exec() == QDialog.DialogCode.Accepted:
            features = list(record.features)
            features[row] = dialog.feature()
            self._replace_current_record(
                features=tuple(features),
                operation="edit_feature",
                parameters={"index": row},
            )

    def delete_selected_features(self) -> None:
        record = self.current_record
        if record is None:
            return
        rows = sorted({index.row() for index in self.feature_table.selectedIndexes()}, reverse=True)
        if not rows:
            return
        features = list(record.features)
        for row in rows:
            if 0 <= row < len(features):
                del features[row]
        self._replace_current_record(
            features=tuple(features),
            operation="delete_feature",
            parameters={"count": len(rows)},
        )

    def add_primer(self) -> None:
        record = self.current_record
        if record is None:
            return
        dialog = PrimerDialog(self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            primer = dialog.primer()
            self._replace_current_record(
                primers=record.primers + (primer,),
                operation="add_primer",
                parameters={"name": primer.name},
            )

    def undo(self) -> None:
        """Restore the previous record snapshot for the current selection."""

        record = self.current_record
        if record is None or not self._undo_stack:
            return
        previous = self._undo_stack.pop()
        self._redo_stack.append(record)
        self._set_current_record(previous, push_undo=False)

    def redo(self) -> None:
        """Reapply the next record snapshot for the current selection."""

        record = self.current_record
        if record is None or not self._redo_stack:
            return
        next_record = self._redo_stack.pop()
        self._undo_stack.append(record)
        self._set_current_record(next_record, push_undo=False)

    def run_feature_detection(self) -> None:
        """Detect common features and let the user review proposals."""

        record = self.current_record
        if record is None:
            return
        try:
            proposals = detect_features(record)
        except Exception as error:  # pragma: no cover - dialog path.
            self._show_error("Feature Detection Failed", str(error))
            return
        if not proposals:
            QMessageBox.information(self, "Detect Features", "No feature matches were found.")
            return
        dialog = FeatureDetectionDialog(self, proposals)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self._apply_feature_detections(dialog.selected_detections())

    def _apply_feature_detections(self, detections: tuple[FeatureDetection, ...]) -> None:
        record = self.current_record
        if record is None or not detections:
            return
        updated = apply_feature_annotations(record, detections)
        if updated is record:
            return
        self._set_current_record(updated)

    def _record_selection_changed(self, row: int) -> None:
        if row < 0:
            return
        mapped_index = self._project_list_indices[row] if row < len(self._project_list_indices) else row
        if mapped_index == self._current_index:
            return
        self._current_index = mapped_index
        self._last_restriction_analysis = None
        self._last_primer_bindings = {}
        self._last_digest_fragments = None
        self._undo_stack.clear()
        self._redo_stack.clear()
        self._refresh_all()
        self._update_undo_redo_actions()

    def _populate_project_list(self, indices: Iterable[int] | None = None) -> None:
        self._project_list_indices = list(range(len(self._records))) if indices is None else list(indices)
        self.project_list.clear()
        for index in self._project_list_indices:
            record = self._records[index]
            self.project_list.addItem(f"{record.id}  ({record.length} bp, {record.topology.value})")
        if self._project_db is None:
            self.project_status_label.setText(f"{len(self._records)} records loaded.")
        else:
            self.project_status_label.setText(f"{len(self._project_list_indices)} of {len(self._records)} project records shown.")

    def search_project_panel(self) -> None:
        """Filter the project panel using text and sequence search."""

        query = self.project_search_edit.text().strip()
        if not query:
            self.clear_project_search()
            return
        matched_indices: set[int] = set()
        if self._project_db is not None and self._project_record_pks:
            pk_to_index = {pk: index for index, pk in enumerate(self._project_record_pks)}
            for match in self._project_db.search_text(query):
                if match.record_pk in pk_to_index:
                    matched_indices.add(pk_to_index[match.record_pk])
            for match in self._project_db.search_sequence(query):
                if match.record_pk in pk_to_index:
                    matched_indices.add(pk_to_index[match.record_pk])
        else:
            lowered = query.lower()
            for index, record in enumerate(self._records):
                text_values = [
                    record.id,
                    record.name or "",
                    record.description or "",
                    *(feature.name or feature.type for feature in record.features),
                    *(primer.name for primer in record.primers),
                ]
                if any(lowered in value.lower() for value in text_values) or query.upper() in record.sequence:
                    matched_indices.add(index)
        self._populate_project_list(sorted(matched_indices))
        if matched_indices:
            self.project_list.setCurrentRow(0)

    def clear_project_search(self) -> None:
        """Clear project-panel filtering."""

        self.project_search_edit.clear()
        self._populate_project_list()
        if self._records:
            visible_row = self._project_list_indices.index(self._current_index) if self._current_index in self._project_list_indices else 0
            self.project_list.setCurrentRow(visible_row)

    def _set_current_record(self, updated: SequenceRecord, *, push_undo: bool = True) -> None:
        record = self.current_record
        if record is None or self._current_index < 0:
            return
        if push_undo:
            self._undo_stack.append(record)
            self._redo_stack.clear()
        if self._project_db is not None and self._current_index < len(self._project_record_pks):
            self._project_db.update_record(self._project_record_pks[self._current_index], updated)
        self._records[self._current_index] = updated
        for visible_row, record_index in enumerate(self._project_list_indices):
            if record_index == self._current_index:
                self.project_list.item(visible_row).setText(
                    f"{updated.id}  ({updated.length} bp, {updated.topology.value})"
                )
                break
        self._refresh_all()
        self._update_undo_redo_actions()

    def _update_undo_redo_actions(self) -> None:
        if hasattr(self, "undo_action"):
            self.undo_action.setEnabled(bool(self._undo_stack))
        if hasattr(self, "redo_action"):
            self.redo_action.setEnabled(bool(self._redo_stack))

    def _replace_current_record(
        self,
        *,
        features: tuple[Feature, ...] | None = None,
        primers: tuple[Primer, ...] | None = None,
        operation: str,
        parameters: dict[str, Any],
    ) -> None:
        record = self.current_record
        if record is None:
            return
        event = ProvenanceEvent(
            operation=operation,
            parameters=parameters,
            input_record_id=record.id,
            input_record_ids=(record.id,),
            input_version_id=record.version_id,
            input_version_ids=(record.version_id,),
            output_record_id=record.id,
            affected_ranges=((0, record.length),),
            description=f"{operation.replace('_', ' ')} on {record.id}",
        )
        updated = record._copy(
            features=None if features is None else features,
            primers=None if primers is None else primers,
            event=event,
        )
        self._set_current_record(updated)

    def _refresh_all(self) -> None:
        self._refresh_map()
        self._refresh_sequence()
        self._refresh_features()
        self._refresh_enzymes()
        self._refresh_primers()
        self._refresh_gel()
        self._refresh_history()

    def _refresh_map(self) -> None:
        record = self.current_record
        if record is None:
            return
        self.map_widget.set_model(render_plasmid_map(record, overlays=self._map_overlays(record)))

    def _map_overlays(self, record: SequenceRecord):
        return prepare_map_overlays(
            record,
            restriction_analysis=self._last_restriction_analysis,
            primer_bindings_by_primer_index=self._last_primer_bindings,
        )

    def _refresh_sequence(self) -> None:
        record = self.current_record
        if record is None:
            self.sequence_view.clear()
            self.translation_view.clear()
            return
        plain_text = _wrapped_sequence_text(record.sequence)
        self.sequence_view.setPlainText(plain_text)
        self._highlight_features(record)
        translations: list[str] = []
        for feature in record.features:
            if feature.type != "CDS":
                continue
            try:
                translations.append(f"{feature.name or feature.type}: {translate_feature(record, feature)}")
            except Exception as error:
                translations.append(f"{feature.name or feature.type}: translation failed ({error})")
        self.translation_view.setPlainText("\n".join(translations) or "No CDS features.")
        self._update_selection_status()

    def _highlight_features(self, record: SequenceRecord) -> None:
        selections: list[Any] = []
        for feature in record.features:
            color = FEATURE_COLORS.get(feature.type, FEATURE_COLORS["misc_feature"])
            for segment in feature.segments:
                for start, end in _line_sliced_ranges(segment.start, segment.end):
                    selection = QTextEdit.ExtraSelection()
                    selection.format = QTextCharFormat()
                    selection.format.setBackground(color)
                    cursor = self.sequence_view.textCursor()
                    cursor.setPosition(_sequence_coord_to_document_pos(start))
                    cursor.setPosition(_sequence_coord_to_document_pos(end), QTextCursor.MoveMode.KeepAnchor)
                    selection.cursor = cursor
                    selections.append(selection)
        self.sequence_view.setExtraSelections(selections)

    def _refresh_features(self) -> None:
        record = self.current_record
        self.feature_table.setRowCount(0)
        if record is None:
            return
        self.feature_table.setRowCount(len(record.features))
        for row, feature in enumerate(record.features):
            self.feature_table.setItem(row, 0, QTableWidgetItem(feature.name or ""))
            self.feature_table.setItem(row, 1, QTableWidgetItem(feature.type))
            self.feature_table.setItem(row, 2, QTableWidgetItem(str(feature.strand)))
            self.feature_table.setItem(row, 3, QTableWidgetItem(str(feature.start)))
            self.feature_table.setItem(row, 4, QTableWidgetItem(str(feature.end)))
            self.feature_table.setItem(row, 5, QTableWidgetItem(json.dumps(dict(feature.qualifiers), sort_keys=True)))

    def _refresh_enzymes(self) -> None:
        self.enzyme_table.setRowCount(0)
        self.enzyme_summary.setText("No enzyme search has been run.")
        if self._last_restriction_analysis is not None:
            self._populate_enzyme_table(self._last_restriction_analysis)

    def _populate_enzyme_table(self, analysis: RestrictionAnalysis) -> None:
        self.enzyme_table.setRowCount(len(analysis.sites))
        for row, site in enumerate(analysis.sites):
            values = (
                site.enzyme_name,
                str(site.start),
                str(site.end),
                str(site.strand),
                str(site.top_cut),
                str(site.bottom_cut),
                site.recognition_sequence,
            )
            for column, value in enumerate(values):
                self.enzyme_table.setItem(row, column, QTableWidgetItem(value))
        self.enzyme_summary.setText(
            "Single cutters: "
            + (", ".join(analysis.single_cutters) or "none")
            + "    Noncutters: "
            + (", ".join(analysis.non_cutters) or "none")
        )

    def _refresh_primers(self) -> None:
        record = self.current_record
        self.primer_table.setRowCount(0)
        self._last_primer_bindings = {}
        if record is None:
            return
        self.primer_table.setRowCount(len(record.primers))
        for row, primer in enumerate(record.primers):
            metrics = primer_metrics(primer.sequence)
            try:
                bindings = find_primer_bindings(record, primer.sequence)
            except Exception:
                bindings = ()
            self._last_primer_bindings[row] = bindings
            binding_text = "; ".join(
                f"{binding.start}-{binding.end} ({'+' if binding.strand == 1 else '-'})"
                for binding in bindings
            )
            values = (
                primer.name,
                primer.sequence,
                str(metrics.length),
                f"{metrics.gc_percent:.1f}",
                f"{metrics.wallace_tm:.1f}",
                str(len(bindings)),
                binding_text or "none",
            )
            for column, value in enumerate(values):
                self.primer_table.setItem(row, column, QTableWidgetItem(value))
        self._refresh_map()

    def _refresh_gel(self) -> None:
        record = self.current_record
        if record is None or record.length == 0:
            return
        fragments: tuple[Any, ...]
        lane_name = "Uncut"
        if self._last_digest_fragments:
            fragments = self._last_digest_fragments
            lane_name = "Digest"
        else:
            fragments = (GelInputFragment(size_bp=record.length, name="uncut"),)
        self.gel_widget.set_model(
            simulate_gel(
                (GelLaneInput(name=lane_name, fragments=fragments),),
                ladder="1kb_dna_ladder",
                amount_mode="equal_mass",
            )
        )

    def _refresh_history(self) -> None:
        record = self.current_record
        self.history_table.setRowCount(0)
        self.history_graph_widget.set_graph(build_history_graph(tuple(self._records)))
        if record is None:
            return
        self.history_table.setRowCount(len(record.history))
        for row, event in enumerate(record.history):
            values = (
                event.timestamp or "",
                event.operation,
                ", ".join(event.input_record_ids) or (event.input_record_id or ""),
                event.output_record_id or "",
                "; ".join(f"{start}-{end}" for start, end in event.affected_ranges),
                event.description or json.dumps(dict(event.parameters), sort_keys=True),
            )
            for column, value in enumerate(values):
                self.history_table.setItem(row, column, QTableWidgetItem(value))

    def _update_selection_status(self) -> None:
        selected = self.sequence_view.textCursor().selectedText()
        bases = [base for base in selected if base.upper() in {"A", "C", "G", "T", "U", "N"}]
        if bases:
            self.selection_label.setText(f"Selection: {len(bases)} bases")
        else:
            self.selection_label.setText("Selection: none")

    def _show_error(self, title: str, message: str) -> None:
        QMessageBox.critical(self, title, message)


class FeatureDialog(QDialog):
    """Small feature editor dialog for simple single-segment features."""

    def __init__(self, parent: QWidget, *, sequence_length: int, feature: Feature | None = None) -> None:
        super().__init__(parent)
        if feature is not None and not _feature_editable_in_simple_dialog(feature):
            msg = "compound features require advanced editing and are not editable in this dialog yet"
            raise ValueError(msg)
        self.setWindowTitle("Feature")
        layout = QFormLayout(self)
        self.name_edit = QLineEdit(feature.name if feature and feature.name else "")
        self.type_edit = QLineEdit(feature.type if feature else "misc_feature")
        self.strand_combo = QComboBox()
        self.strand_combo.addItems(["1", "-1", "0"])
        if feature:
            self.strand_combo.setCurrentText(str(feature.strand))
        self.start_spin = QSpinBox()
        self.start_spin.setRange(0, sequence_length)
        self.end_spin = QSpinBox()
        self.end_spin.setRange(0, sequence_length)
        if feature:
            self.start_spin.setValue(feature.start)
            self.end_spin.setValue(feature.end)
        else:
            self.end_spin.setValue(min(sequence_length, 1))
        self.qualifiers_edit = QLineEdit(
            json.dumps(dict(feature.qualifiers), sort_keys=True) if feature else "{}"
        )
        layout.addRow("Name", self.name_edit)
        layout.addRow("Type", self.type_edit)
        layout.addRow("Strand", self.strand_combo)
        layout.addRow("Start", self.start_spin)
        layout.addRow("End", self.end_spin)
        layout.addRow("Qualifiers JSON", self.qualifiers_edit)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addRow(buttons)

    def feature(self) -> Feature:
        qualifiers = json.loads(self.qualifiers_edit.text() or "{}")
        return Feature(
            type=self.type_edit.text(),
            name=self.name_edit.text() or None,
            start=self.start_spin.value(),
            end=self.end_spin.value(),
            strand=int(self.strand_combo.currentText()),
            qualifiers=qualifiers,
        )


def _feature_editable_in_simple_dialog(feature: Feature) -> bool:
    return len(feature.segments) == 1


class PrimerDialog(QDialog):
    """Small primer add dialog."""

    def __init__(self, parent: QWidget) -> None:
        super().__init__(parent)
        self.setWindowTitle("Primer")
        layout = QFormLayout(self)
        self.name_edit = QLineEdit("primer")
        self.sequence_edit = QLineEdit()
        layout.addRow("Name", self.name_edit)
        layout.addRow("Sequence", self.sequence_edit)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addRow(buttons)

    def primer(self) -> Primer:
        return Primer(name=self.name_edit.text(), sequence=self.sequence_edit.text())


class FeatureDetectionDialog(QDialog):
    """Review proposed automatic feature annotations before applying them."""

    def __init__(self, parent: QWidget, detections: tuple[FeatureDetection, ...]) -> None:
        super().__init__(parent)
        self._detections = tuple(detections)
        self.setWindowTitle("Detect Features")
        self.resize(900, 420)
        layout = QVBoxLayout(self)
        self.table = QTableWidget(len(self._detections), 7)
        self.table.setHorizontalHeaderLabels(
            ["Apply", "Name", "Type", "Range", "Strand", "Identity", "Notes"]
        )
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        for row, detection in enumerate(self._detections):
            apply_item = QTableWidgetItem()
            apply_item.setFlags(
                Qt.ItemFlag.ItemIsEnabled
                | Qt.ItemFlag.ItemIsUserCheckable
                | Qt.ItemFlag.ItemIsSelectable
            )
            apply_item.setCheckState(Qt.CheckState.Checked)
            self.table.setItem(row, 0, apply_item)
            values = (
                detection.name,
                detection.type,
                _detection_range_text(detection),
                {1: "+", -1: "-", 0: "unknown"}[detection.strand],
                f"{detection.identity * 100:.1f}%",
                detection.entry.notes,
            )
            for column, value in enumerate(values, start=1):
                item = QTableWidgetItem(value)
                item.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
                self.table.setItem(row, column, item)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(self.table)
        layout.addWidget(buttons)

    def selected_detections(self) -> tuple[FeatureDetection, ...]:
        """Return proposals the user left checked."""

        selected: list[FeatureDetection] = []
        for row, detection in enumerate(self._detections):
            item = self.table.item(row, 0)
            if item is not None and item.checkState() == Qt.CheckState.Checked:
                selected.append(detection)
        return tuple(selected)

def _detection_range_text(detection: FeatureDetection) -> str:
    return ", ".join(f"{segment.start}-{segment.end}" for segment in detection.segments)


def _wrapped_sequence_text(sequence: str, *, width: int = 60) -> str:
    lines = []
    for start in range(0, len(sequence), width):
        lines.append(f"{start + 1:>8}  {sequence[start:start + width]}")
    return "\n".join(lines)


def _sequence_coord_to_document_pos(coord: int, *, width: int = 60, prefix: int = 10) -> int:
    line = coord // width
    column = coord % width
    return line * (prefix + width + 1) + prefix + column


def _line_sliced_ranges(start: int, end: int, *, width: int = 60) -> list[tuple[int, int]]:
    ranges: list[tuple[int, int]] = []
    cursor = start
    while cursor < end:
        line_end = min(end, ((cursor // width) + 1) * width)
        ranges.append((cursor, line_end))
        cursor = line_end
    return ranges


def main(argv: list[str] | None = None) -> int:
    """Open the PlasmidLab desktop application."""

    arguments = list(sys.argv[1:] if argv is None else argv)
    if "--help" in arguments or "-h" in arguments:
        _write_stdout(_gui_cli_help())
        return 0
    if "--version" in arguments:
        import plasmidlab

        _write_stdout(f"PlasmidLab {plasmidlab.__version__}\n")
        return 0
    if "--smoke-json" in arguments:
        index = arguments.index("--smoke-json")
        try:
            output_path = Path(arguments[index + 1])
        except IndexError as error:
            msg = "--smoke-json requires an output path"
            raise SystemExit(msg) from error
        from importlib import resources

        import plasmidlab
        from plasmidlab.core.feature_annotation import load_feature_library

        feature_resource_root = resources.files("plasmidlab.data.features")
        resource_entries = [
            resource for resource in feature_resource_root.iterdir() if resource.name.endswith(".json")
        ]
        library = load_feature_library()
        package_resources_available = bool(resource_entries)
        feature_library_available = bool(library)
        output_path.write_text(
            json.dumps(
                {
                    "application": "PlasmidLab",
                    "cli": True,
                    "feature_library": feature_library_available,
                    "version": plasmidlab.__version__,
                    "package_resources": package_resources_available,
                    "feature_library_entries": len(library),
                    "feature_resource_entries": len(resource_entries),
                    "ok": package_resources_available and feature_library_available,
                    "python_version_ok": sys.version_info[:2] >= (3, 12),
                    "qt_import": QApplication is not None,
                },
                sort_keys=True,
            ),
            encoding="utf-8",
        )
        return 0

    if QApplication is None:
        msg = "PySide6 is required to run the PlasmidLab GUI"
        raise ImportError(msg)
    app = QApplication.instance() or QApplication(sys.argv)
    window = PlasmidLabMainWindow()
    window.show()
    return app.exec()


def _gui_cli_help() -> str:
    return (
        "usage: PlasmidLab [--help] [--version] [--smoke-json PATH]\n\n"
        "Open the PlasmidLab desktop GUI.\n\n"
        "options:\n"
        "  --help             show this help and exit\n"
        "  --version          print the PlasmidLab version and exit\n"
        "  --smoke-json PATH  write a headless frozen-app smoke report and exit\n"
    )


def _write_stdout(text: str) -> None:
    stream = getattr(sys, "stdout", None)
    if stream is not None:
        stream.write(text)


if __name__ == "__main__":
    raise SystemExit(main())
