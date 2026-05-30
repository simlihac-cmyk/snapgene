"""SQLite-backed project database, search, and batch operations."""

from __future__ import annotations

import json
import re
import sqlite3
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from plasmidlab.core import (
    Feature,
    FeatureSegment,
    MoleculeTopology,
    MoleculeType,
    Primer,
    ProvenanceEvent,
    SequenceRecord,
    analyze_restriction_sites,
    apply_feature_annotations,
    detect_features,
)
from plasmidlab.io.fasta import read_fasta, write_fasta
from plasmidlab.io.genbank import read_genbank, write_genbank
from plasmidlab.render import render_plasmid_map, write_svg


DNA_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
GENBANK_SUFFIXES = {".gb", ".gbk", ".genbank"}
FASTA_SUFFIXES = {".fa", ".fasta", ".fna"}


@dataclass(frozen=True, slots=True)
class ProjectTextMatch:
    """A full-text search hit in a project collection."""

    record_pk: int
    record_id: str
    field: str
    text: str


@dataclass(frozen=True, slots=True)
class ProjectSequenceMatch:
    """A sequence-query hit in a project collection."""

    record_pk: int
    record_id: str
    start: int
    end: int
    strand: int
    mismatches: int
    matched_sequence: str
    wraps_origin: bool = False


class ProjectDatabase:
    """SQLite-backed PlasmidLab project collection."""

    def __init__(self, path: str | Path = ":memory:") -> None:
        self.path = str(path)
        self.connection = sqlite3.connect(self.path)
        self.connection.row_factory = sqlite3.Row
        self.connection.execute("PRAGMA foreign_keys = ON")
        self._fts_enabled = False
        self._initialize_schema()

    def close(self) -> None:
        """Close the SQLite connection."""

        self.connection.close()

    def add_record(
        self,
        record: SequenceRecord,
        *,
        tags: Iterable[str] = (),
        note: str | None = None,
        replace: bool = False,
        source_path: str | Path | None = None,
    ) -> int:
        """Add a sequence record and return its database primary key."""

        existing = self.record_pk(record.id)
        if existing is not None:
            if not replace:
                msg = f"record already exists in project: {record.id}"
                raise ValueError(msg)
            self.update_record(existing, record, tags=tags, notes=(note,) if note else None, source_path=source_path)
            return existing

        with self.connection:
            cursor = self.connection.execute(
                """
                INSERT INTO sequences
                    (record_id, version_id, name, description, sequence, molecule_type,
                     topology, length, source_path)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    record.id,
                    record.version_id,
                    record.name,
                    record.description,
                    record.sequence,
                    record.molecule_type.value,
                    record.topology.value,
                    record.length,
                    None if source_path is None else str(source_path),
                ),
            )
            record_pk = int(cursor.lastrowid)
            self._insert_record_children(record_pk, record)
            self.set_tags(record_pk, tags)
            if note:
                self._insert_note(record_pk, note)
            self._index_record(record_pk)
        return record_pk

    def update_record(
        self,
        record: int | str,
        updated: SequenceRecord,
        *,
        tags: Iterable[str] | None = None,
        notes: Iterable[str] | None = None,
        source_path: str | Path | None = None,
    ) -> int:
        """Replace a stored record while preserving its database primary key."""

        record_pk = self._resolve_record_pk(record)
        with self.connection:
            self.connection.execute(
                """
                UPDATE sequences
                SET record_id = ?, name = ?, description = ?, sequence = ?, molecule_type = ?,
                    topology = ?, length = ?, version_id = ?, source_path = COALESCE(?, source_path)
                WHERE id = ?
                """,
                (
                    updated.id,
                    updated.name,
                    updated.description,
                    updated.sequence,
                    updated.molecule_type.value,
                    updated.topology.value,
                    updated.length,
                    updated.version_id,
                    None if source_path is None else str(source_path),
                    record_pk,
                ),
            )
            self._delete_record_children(record_pk)
            self._insert_record_children(record_pk, updated)
            if tags is not None:
                self.set_tags(record_pk, tags)
            if notes is not None:
                self.connection.execute("DELETE FROM notes WHERE sequence_id = ?", (record_pk,))
                for note in notes:
                    self._insert_note(record_pk, note)
            self._delete_search_index(record_pk)
            self._index_record(record_pk)
        return record_pk

    def record_pk(self, record_id: str) -> int | None:
        """Return the database id for a record id, if present."""

        row = self.connection.execute(
            "SELECT id FROM sequences WHERE record_id = ?",
            (record_id,),
        ).fetchone()
        return None if row is None else int(row["id"])

    def get_record(self, record: int | str) -> SequenceRecord:
        """Load a record by database id or record id."""

        record_pk = self._resolve_record_pk(record)
        row = self.connection.execute("SELECT * FROM sequences WHERE id = ?", (record_pk,)).fetchone()
        if row is None:
            msg = f"unknown project record: {record}"
            raise KeyError(msg)
        return SequenceRecord(
            id=row["record_id"],
            name=row["name"],
            description=row["description"],
            sequence=row["sequence"],
            molecule_type=row["molecule_type"],
            topology=row["topology"],
            version_id=row["version_id"],
            features=self._features_for_record(record_pk),
            primers=self._primers_for_record(record_pk),
            history=self._events_for_record(record_pk),
        )

    def all_records(self) -> tuple[SequenceRecord, ...]:
        """Return all records in insertion order."""

        return tuple(self.get_record(int(row["id"])) for row in self.connection.execute("SELECT id FROM sequences ORDER BY id"))

    def record_pks(self) -> tuple[int, ...]:
        """Return database primary keys in project-list order."""

        return tuple(int(row["id"]) for row in self.connection.execute("SELECT id FROM sequences ORDER BY id"))

    def add_note(self, record: int | str, text: str) -> int:
        """Add a note to a record."""

        record_pk = self._resolve_record_pk(record)
        cursor = self._insert_note(record_pk, text)
        self._delete_search_index(record_pk)
        self._index_record(record_pk)
        return int(cursor.lastrowid)

    def notes(self, record: int | str) -> tuple[str, ...]:
        """Return notes for a record."""

        record_pk = self._resolve_record_pk(record)
        return tuple(
            str(row["text"])
            for row in self.connection.execute("SELECT text FROM notes WHERE sequence_id = ? ORDER BY id", (record_pk,))
        )

    def set_tags(self, record: int | str, tags: Iterable[str]) -> None:
        """Replace tags for a record."""

        record_pk = self._resolve_record_pk(record) if not isinstance(record, int) else record
        self.connection.execute("DELETE FROM sequence_tags WHERE sequence_id = ?", (record_pk,))
        for tag in tags:
            normalized = str(tag).strip()
            if not normalized:
                continue
            cursor = self.connection.execute(
                "INSERT OR IGNORE INTO tags (name) VALUES (?)",
                (normalized,),
            )
            tag_row = self.connection.execute("SELECT id FROM tags WHERE name = ?", (normalized,)).fetchone()
            tag_id = int(tag_row["id"] if tag_row is not None else cursor.lastrowid)
            self.connection.execute(
                "INSERT OR IGNORE INTO sequence_tags (sequence_id, tag_id) VALUES (?, ?)",
                (record_pk, tag_id),
            )

    def tags(self, record: int | str) -> tuple[str, ...]:
        """Return tags for a record."""

        record_pk = self._resolve_record_pk(record)
        return tuple(
            str(row["name"])
            for row in self.connection.execute(
                """
                SELECT tags.name
                FROM tags
                JOIN sequence_tags ON sequence_tags.tag_id = tags.id
                WHERE sequence_tags.sequence_id = ?
                ORDER BY tags.name
                """,
                (record_pk,),
            )
        )

    def search_text(self, query: str) -> tuple[ProjectTextMatch, ...]:
        """Full-text search sequence names, feature names, primer names, and notes."""

        query = query.strip()
        if not query:
            return ()
        if self._fts_enabled:
            rows = self.connection.execute(
                """
                SELECT search_index.record_pk, sequences.record_id, search_index.field, search_index.content
                FROM search_index
                JOIN sequences ON sequences.id = search_index.record_pk
                WHERE search_index MATCH ?
                ORDER BY sequences.record_id, search_index.field
                """,
                (_fts_query(query),),
            ).fetchall()
        else:
            rows = self.connection.execute(
                """
                SELECT search_index.record_pk, sequences.record_id, search_index.field, search_index.content
                FROM search_index
                JOIN sequences ON sequences.id = search_index.record_pk
                WHERE lower(search_index.content) LIKE ?
                ORDER BY sequences.record_id, search_index.field
                """,
                (f"%{query.lower()}%",),
            ).fetchall()
        return tuple(
            ProjectTextMatch(
                record_pk=int(row["record_pk"]),
                record_id=str(row["record_id"]),
                field=str(row["field"]),
                text=str(row["content"]),
            )
            for row in rows
        )

    def search_sequence(
        self,
        query: str,
        *,
        record_ids: Iterable[int | str] | None = None,
        mismatches: int = 0,
        include_reverse_complement: bool = True,
    ) -> tuple[ProjectSequenceMatch, ...]:
        """Search records for a DNA or protein subsequence."""

        normalized_query = "".join(query.upper().split()).replace("U", "T")
        if not normalized_query:
            return ()
        if mismatches < 0:
            msg = "mismatches must be non-negative"
            raise ValueError(msg)
        records = self._selected_record_rows(record_ids)
        matches: list[ProjectSequenceMatch] = []
        for row in records:
            record_pk = int(row["id"])
            record_id = str(row["record_id"])
            sequence = str(row["sequence"]).upper()
            topology = MoleculeTopology(row["topology"])
            molecule_type = MoleculeType(row["molecule_type"])
            patterns = [(normalized_query, 1)]
            if include_reverse_complement and molecule_type is MoleculeType.DNA and _looks_like_dna(normalized_query):
                reverse = _reverse_complement(normalized_query)
                if reverse != normalized_query:
                    patterns.append((reverse, -1))
            seen: set[tuple[int, int, int]] = set()
            for pattern, strand in patterns:
                for start, end, wraps, matched, mismatch_count in _find_sequence_matches(
                    sequence,
                    pattern,
                    topology,
                    mismatches,
                ):
                    key = (start, end, strand)
                    if key in seen:
                        continue
                    seen.add(key)
                    matches.append(
                        ProjectSequenceMatch(
                            record_pk=record_pk,
                            record_id=record_id,
                            start=start,
                            end=end,
                            strand=strand,
                            mismatches=mismatch_count,
                            matched_sequence=matched,
                            wraps_origin=wraps,
                        )
                    )
        return tuple(sorted(matches, key=lambda item: (item.record_id, item.start, item.strand, item.mismatches)))

    def import_folder(
        self,
        folder: str | Path,
        *,
        recursive: bool = False,
        replace: bool = False,
    ) -> tuple[int, ...]:
        """Import all FASTA and GenBank records from a folder."""

        root = Path(folder)
        pattern = "**/*" if recursive else "*"
        paths = sorted(
            {
                path
                for path in root.glob(pattern)
                if path.is_file() and path.suffix.lower() in GENBANK_SUFFIXES | FASTA_SUFFIXES
            },
            key=lambda path: str(path),
        )
        imported: list[int] = []
        for path in paths:
            if path.suffix.lower() in GENBANK_SUFFIXES:
                records = read_genbank(path)
            else:
                records = read_fasta(path)
            for record in records:
                imported.append(self.add_record(record, replace=replace, source_path=path))
        return tuple(imported)

    def export_records(
        self,
        records: Iterable[int | str],
        output_dir: str | Path,
        *,
        file_format: str = "genbank",
    ) -> tuple[Path, ...]:
        """Export selected records to FASTA or GenBank files."""

        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)
        exported: list[Path] = []
        for record_key in records:
            record = self.get_record(record_key)
            if file_format.lower() in {"fa", "fasta"}:
                path = output / f"{_safe_filename(record.id)}.fasta"
                write_fasta(record, path)
            elif file_format.lower() in {"gb", "gbk", "genbank"}:
                path = output / f"{_safe_filename(record.id)}.gb"
                write_genbank(record, path)
            else:
                msg = f"unsupported export format: {file_format}"
                raise ValueError(msg)
            exported.append(path)
        return tuple(exported)

    def batch_detect_features(
        self,
        records: Iterable[int | str],
        *,
        library_paths: str | Path | Iterable[str | Path] | None = None,
    ) -> dict[int, int]:
        """Detect and apply common features to selected records."""

        applied: dict[int, int] = {}
        for record_key in records:
            record_pk = self._resolve_record_pk(record_key)
            record = self.get_record(record_pk)
            detections = detect_features(record, library_paths=library_paths)
            updated = apply_feature_annotations(record, detections)
            if updated is not record:
                self.update_record(record_pk, updated)
            applied[record_pk] = len(updated.features) - len(record.features)
        return applied

    def batch_restriction_analysis(
        self,
        records: Iterable[int | str],
        enzymes: Any,
    ) -> dict[int, Any]:
        """Run restriction analysis on selected records."""

        return {
            self._resolve_record_pk(record_key): analyze_restriction_sites(self.get_record(record_key), enzymes)
            for record_key in records
        }

    def batch_export_map_svg(
        self,
        records: Iterable[int | str],
        output_dir: str | Path,
    ) -> tuple[Path, ...]:
        """Export plasmid/sequence maps for selected records as SVG."""

        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)
        exported: list[Path] = []
        for record_key in records:
            record = self.get_record(record_key)
            path = output / f"{_safe_filename(record.id)}.svg"
            write_svg(render_plasmid_map(record), path)
            exported.append(path)
        return tuple(exported)

    def _initialize_schema(self) -> None:
        with self.connection:
            self.connection.executescript(
                """
                CREATE TABLE IF NOT EXISTS sequences (
                    id INTEGER PRIMARY KEY,
                    record_id TEXT NOT NULL UNIQUE,
                    version_id TEXT,
                    name TEXT,
                    description TEXT,
                    sequence TEXT NOT NULL,
                    molecule_type TEXT NOT NULL,
                    topology TEXT NOT NULL,
                    length INTEGER NOT NULL,
                    source_path TEXT,
                    created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
                    updated_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
                );
                CREATE TABLE IF NOT EXISTS features (
                    id INTEGER PRIMARY KEY,
                    sequence_id INTEGER NOT NULL REFERENCES sequences(id) ON DELETE CASCADE,
                    name TEXT,
                    type TEXT NOT NULL,
                    start INTEGER NOT NULL,
                    end INTEGER NOT NULL,
                    strand INTEGER NOT NULL,
                    segments_json TEXT NOT NULL,
                    qualifiers_json TEXT NOT NULL
                );
                CREATE TABLE IF NOT EXISTS primers (
                    id INTEGER PRIMARY KEY,
                    sequence_id INTEGER NOT NULL REFERENCES sequences(id) ON DELETE CASCADE,
                    name TEXT NOT NULL,
                    sequence TEXT NOT NULL,
                    start INTEGER,
                    end INTEGER,
                    strand INTEGER NOT NULL,
                    target_id TEXT
                );
                CREATE TABLE IF NOT EXISTS provenance_events (
                    id INTEGER PRIMARY KEY,
                    sequence_id INTEGER NOT NULL REFERENCES sequences(id) ON DELETE CASCADE,
                    event_order INTEGER NOT NULL,
                    operation TEXT NOT NULL,
                    timestamp TEXT,
                    input_record_id TEXT,
                    input_record_ids_json TEXT NOT NULL,
                    output_record_id TEXT,
                    input_version_id TEXT,
                    input_version_ids_json TEXT NOT NULL DEFAULT '[]',
                    output_version_id TEXT,
                    parameters_json TEXT NOT NULL,
                    affected_ranges_json TEXT NOT NULL,
                    description TEXT NOT NULL
                );
                CREATE TABLE IF NOT EXISTS tags (
                    id INTEGER PRIMARY KEY,
                    name TEXT NOT NULL UNIQUE
                );
                CREATE TABLE IF NOT EXISTS sequence_tags (
                    sequence_id INTEGER NOT NULL REFERENCES sequences(id) ON DELETE CASCADE,
                    tag_id INTEGER NOT NULL REFERENCES tags(id) ON DELETE CASCADE,
                    PRIMARY KEY (sequence_id, tag_id)
                );
                CREATE TABLE IF NOT EXISTS notes (
                    id INTEGER PRIMARY KEY,
                    sequence_id INTEGER NOT NULL REFERENCES sequences(id) ON DELETE CASCADE,
                    text TEXT NOT NULL,
                    created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
                );
                """
            )
            self._ensure_schema_migrations()
            try:
                self.connection.execute(
                    "CREATE VIRTUAL TABLE IF NOT EXISTS search_index USING fts5(record_pk UNINDEXED, field UNINDEXED, content)"
                )
                self._fts_enabled = True
            except sqlite3.OperationalError:
                self.connection.execute(
                    """
                    CREATE TABLE IF NOT EXISTS search_index (
                        record_pk INTEGER NOT NULL,
                        field TEXT NOT NULL,
                        content TEXT NOT NULL
                    )
                    """
                )
                self._fts_enabled = False

    def _ensure_schema_migrations(self) -> None:
        sequence_columns = _table_columns(self.connection, "sequences")
        if "version_id" not in sequence_columns:
            self.connection.execute("ALTER TABLE sequences ADD COLUMN version_id TEXT")

        event_columns = _table_columns(self.connection, "provenance_events")
        if "input_version_id" not in event_columns:
            self.connection.execute("ALTER TABLE provenance_events ADD COLUMN input_version_id TEXT")
        if "input_version_ids_json" not in event_columns:
            self.connection.execute(
                "ALTER TABLE provenance_events ADD COLUMN input_version_ids_json TEXT NOT NULL DEFAULT '[]'"
            )
        if "output_version_id" not in event_columns:
            self.connection.execute("ALTER TABLE provenance_events ADD COLUMN output_version_id TEXT")

    def _insert_record_children(self, record_pk: int, record: SequenceRecord) -> None:
        for feature in record.features:
            self.connection.execute(
                """
                INSERT INTO features
                    (sequence_id, name, type, start, end, strand, segments_json, qualifiers_json)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    record_pk,
                    feature.name,
                    feature.type,
                    feature.start,
                    feature.end,
                    feature.strand,
                    json.dumps([(segment.start, segment.end, segment.strand) for segment in feature.segments]),
                    json.dumps(dict(feature.qualifiers), sort_keys=True),
                ),
            )
        for primer in record.primers:
            self.connection.execute(
                """
                INSERT INTO primers (sequence_id, name, sequence, start, end, strand, target_id)
                VALUES (?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    record_pk,
                    primer.name,
                    primer.sequence,
                    primer.start,
                    primer.end,
                    primer.strand,
                    primer.target_id,
                ),
            )
        for order, event in enumerate(record.history):
            self.connection.execute(
                """
                INSERT INTO provenance_events
                    (sequence_id, event_order, operation, timestamp, input_record_id,
                     input_record_ids_json, output_record_id, input_version_id,
                     input_version_ids_json, output_version_id, parameters_json,
                     affected_ranges_json, description)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    record_pk,
                    order,
                    event.operation,
                    event.timestamp,
                    event.input_record_id,
                    json.dumps(event.input_record_ids),
                    event.output_record_id,
                    event.input_version_id,
                    json.dumps(event.input_version_ids),
                    event.output_version_id,
                    json.dumps(_json_safe(dict(event.parameters)), sort_keys=True),
                    json.dumps(event.affected_ranges),
                    event.description,
                ),
            )

    def _delete_record_children(self, record_pk: int) -> None:
        for table in ("features", "primers", "provenance_events"):
            self.connection.execute(f"DELETE FROM {table} WHERE sequence_id = ?", (record_pk,))

    def _insert_note(self, record_pk: int, text: str) -> sqlite3.Cursor:
        return self.connection.execute(
            "INSERT INTO notes (sequence_id, text) VALUES (?, ?)",
            (record_pk, text),
        )

    def _features_for_record(self, record_pk: int) -> tuple[Feature, ...]:
        features: list[Feature] = []
        for row in self.connection.execute("SELECT * FROM features WHERE sequence_id = ? ORDER BY id", (record_pk,)):
            segments = tuple(
                FeatureSegment(int(start), int(end), int(strand))
                for start, end, strand in json.loads(row["segments_json"])
            )
            features.append(
                Feature(
                    type=row["type"],
                    name=row["name"],
                    segments=segments,
                    qualifiers=json.loads(row["qualifiers_json"]),
                )
            )
        return tuple(features)

    def _primers_for_record(self, record_pk: int) -> tuple[Primer, ...]:
        return tuple(
            Primer(
                name=row["name"],
                sequence=row["sequence"],
                start=row["start"],
                end=row["end"],
                strand=int(row["strand"]),
                target_id=row["target_id"],
            )
            for row in self.connection.execute("SELECT * FROM primers WHERE sequence_id = ? ORDER BY id", (record_pk,))
        )

    def _events_for_record(self, record_pk: int) -> tuple[ProvenanceEvent, ...]:
        return tuple(
            ProvenanceEvent(
                operation=row["operation"],
                parameters=json.loads(row["parameters_json"]),
                input_record_id=row["input_record_id"],
                timestamp=row["timestamp"],
                input_record_ids=tuple(json.loads(row["input_record_ids_json"])),
                output_record_id=row["output_record_id"],
                input_version_id=row["input_version_id"],
                input_version_ids=tuple(json.loads(row["input_version_ids_json"] or "[]")),
                output_version_id=row["output_version_id"],
                affected_ranges=tuple(tuple(item) for item in json.loads(row["affected_ranges_json"])),
                description=row["description"],
            )
            for row in self.connection.execute(
                "SELECT * FROM provenance_events WHERE sequence_id = ? ORDER BY event_order",
                (record_pk,),
            )
        )

    def _index_record(self, record_pk: int) -> None:
        record = self.get_record(record_pk)
        values = [
            ("sequence", " ".join(item for item in (record.id, record.name or "", record.description or "") if item)),
        ]
        values.extend(("feature", feature.name or feature.type) for feature in record.features)
        values.extend(("primer", primer.name) for primer in record.primers)
        values.extend(("note", note) for note in self.notes(record_pk))
        for field, content in values:
            if content:
                self.connection.execute(
                    "INSERT INTO search_index (record_pk, field, content) VALUES (?, ?, ?)",
                    (record_pk, field, content),
                )

    def _delete_search_index(self, record_pk: int) -> None:
        self.connection.execute("DELETE FROM search_index WHERE record_pk = ?", (record_pk,))

    def _resolve_record_pk(self, record: int | str) -> int:
        if isinstance(record, int):
            row = self.connection.execute("SELECT id FROM sequences WHERE id = ?", (record,)).fetchone()
            if row is not None:
                return record
        else:
            found = self.record_pk(record)
            if found is not None:
                return found
        msg = f"unknown project record: {record}"
        raise KeyError(msg)

    def _selected_record_rows(self, record_ids: Iterable[int | str] | None) -> tuple[sqlite3.Row, ...]:
        if record_ids is None:
            return tuple(self.connection.execute("SELECT * FROM sequences ORDER BY id"))
        rows = []
        for record_key in record_ids:
            record_pk = self._resolve_record_pk(record_key)
            rows.append(self.connection.execute("SELECT * FROM sequences WHERE id = ?", (record_pk,)).fetchone())
        return tuple(row for row in rows if row is not None)


def _find_sequence_matches(
    sequence: str,
    pattern: str,
    topology: MoleculeTopology,
    mismatch_tolerance: int,
) -> Iterable[tuple[int, int, bool, str, int]]:
    if not sequence or not pattern:
        return
    if len(pattern) > len(sequence):
        return
    max_start = len(sequence) if topology is MoleculeTopology.CIRCULAR else len(sequence) - len(pattern) + 1
    for start in range(max_start):
        end_abs = start + len(pattern)
        if end_abs <= len(sequence):
            window = sequence[start:end_abs]
            wraps = False
        elif topology is MoleculeTopology.CIRCULAR:
            window = sequence[start:] + sequence[: end_abs % len(sequence)]
            wraps = True
        else:
            continue
        mismatches = sum(1 for left, right in zip(window, pattern, strict=True) if left != right)
        if mismatches <= mismatch_tolerance:
            yield (start, end_abs % len(sequence) if wraps else end_abs, wraps, window, mismatches)


def _reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(DNA_COMPLEMENT)[::-1]


def _looks_like_dna(sequence: str) -> bool:
    return bool(sequence) and set(sequence.upper()) <= set("ACGTRYSWKMBDHVN")


def _fts_query(query: str) -> str:
    tokens = re.findall(r"[\w.-]+", query)
    return " ".join(f'"{token}"' for token in tokens) if tokens else '""'


def _safe_filename(value: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return sanitized or "record"


def _table_columns(connection: sqlite3.Connection, table: str) -> set[str]:
    return {str(row["name"]) for row in connection.execute(f"PRAGMA table_info({table})")}


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, tuple | list):
        return [_json_safe(item) for item in value]
    return value
