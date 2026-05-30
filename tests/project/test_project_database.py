import sqlite3

from plasmidlab.core import Feature, MoleculeTopology, Primer, SequenceRecord
from plasmidlab.history import build_history_graph
from plasmidlab.io.fasta import write_fasta, read_fasta
from plasmidlab.io.genbank import read_genbank, write_genbank
from plasmidlab.project import ProjectDatabase


def test_project_imports_multiple_genbank_and_fasta_files(tmp_path) -> None:
    gb_record = SequenceRecord(
        id="plasmid_a",
        name="Plasmid_A",
        sequence="A" * 12 + "TAATACGACTCACTATAGGG" + "C" * 12,
        topology=MoleculeTopology.CIRCULAR,
        features=(Feature(type="promoter", start=12, end=32, name="T7 promoter"),),
    )
    fasta_record = SequenceRecord(id="linear_b", name="Linear_B", sequence="GGGAAACCCGGG")
    write_genbank(gb_record, tmp_path / "plasmid_a.gb")
    write_fasta(fasta_record, tmp_path / "linear_b.fasta")

    db = ProjectDatabase(tmp_path / "project.sqlite")
    imported = db.import_folder(tmp_path)

    assert len(imported) == 2
    assert [record.id for record in db.all_records()] == ["linear_b", "plasmid_a"]
    loaded = db.get_record("plasmid_a")
    assert loaded.features[0].name == "T7 promoter"

    primer_record = SequenceRecord(
        id="with_primer",
        sequence="TAATACGACTCA",
        primers=(Primer(name="T7 forward", sequence="TAATACGACTCA"),),
    )
    db.add_record(primer_record)
    assert db.get_record("with_primer").primers[0].name == "T7 forward"


def test_project_full_text_search_by_feature_name_and_notes(tmp_path) -> None:
    record = SequenceRecord(
        id="searchable",
        sequence="ACGTACGT",
        features=(Feature(type="CDS", start=0, end=4, name="AmpR marker"),),
    )
    db = ProjectDatabase(tmp_path / "project.sqlite")
    pk = db.add_record(record, note="favorite construct")

    feature_hits = db.search_text("AmpR")
    note_hits = db.search_text("favorite")

    assert [(hit.record_id, hit.field) for hit in feature_hits] == [("searchable", "feature")]
    assert [(hit.record_pk, hit.field) for hit in note_hits] == [(pk, "note")]


def test_project_sequence_search_exact_reverse_complement_and_mismatch(tmp_path) -> None:
    db = ProjectDatabase(tmp_path / "project.sqlite")
    db.add_record(
        SequenceRecord(
            id="seq",
            sequence="AAACCCGGGTTT",
            topology=MoleculeTopology.CIRCULAR,
        )
    )

    exact = db.search_sequence("CCCGGG")
    reverse = db.search_sequence("AAACCC")
    mismatch = db.search_sequence("CCCGGA", mismatches=1)
    crossing = db.search_sequence("TTTAAA")

    assert [(hit.start, hit.end, hit.strand, hit.mismatches) for hit in exact] == [(3, 9, 1, 0)]
    assert any(hit.strand == -1 for hit in reverse)
    assert mismatch[0].mismatches == 1
    assert crossing[0].wraps_origin
    assert crossing[0].start == 9
    assert crossing[0].end == 3


def test_project_batch_export_round_trip(tmp_path) -> None:
    db = ProjectDatabase(tmp_path / "project.sqlite")
    first = db.add_record(
        SequenceRecord(
            id="first",
            sequence="ATGCATGC",
            features=(Feature(type="misc_feature", start=0, end=4, name="left"),),
        )
    )
    second = db.add_record(SequenceRecord(id="second", sequence="GGGGCCCC"))
    export_dir = tmp_path / "export"

    paths = db.export_records((first, second), export_dir, file_format="genbank")
    loaded = tuple(record for path in paths for record in read_genbank(path))

    assert [path.name for path in paths] == ["first.gb", "second.gb"]
    assert [record.id for record in loaded] == ["first", "second"]
    assert loaded[0].features[0].name == "left"


def test_project_batch_operations(tmp_path) -> None:
    db = ProjectDatabase(tmp_path / "project.sqlite")
    record_pk = db.add_record(
        SequenceRecord(
            id="batch",
            sequence="TAATACGACTCACTATAGGG" + "A" * 20 + "GAATTC",
            topology=MoleculeTopology.CIRCULAR,
        )
    )
    svg_dir = tmp_path / "maps"

    detected = db.batch_detect_features((record_pk,))
    analysis = db.batch_restriction_analysis((record_pk,), "EcoRI")
    svg_paths = db.batch_export_map_svg((record_pk,), svg_dir)

    assert detected[record_pk] >= 1
    assert db.get_record(record_pk).features
    assert len(analysis[record_pk].sites) == 1
    assert svg_paths[0].read_text(encoding="utf-8").startswith('<?xml version="1.0"')


def test_project_fasta_export_round_trip(tmp_path) -> None:
    db = ProjectDatabase(tmp_path / "project.sqlite")
    pk = db.add_record(SequenceRecord(id="fasta_record", sequence="ATGCATGC"))
    paths = db.export_records((pk,), tmp_path / "fasta_export", file_format="fasta")

    loaded = read_fasta(paths[0])

    assert loaded[0].id == "fasta_record"
    assert loaded[0].sequence == "ATGCATGC"


def test_project_db_preserves_version_identity_for_history_graph(tmp_path) -> None:
    edited = SequenceRecord(id="project_versions", sequence="AAAA").insert(1, "GG").delete(0, 1)
    db = ProjectDatabase(tmp_path / "project.sqlite")

    pk = db.add_record(edited)
    reloaded = db.get_record(pk)
    graph = build_history_graph(reloaded)
    record_node_ids = {node.id for node in graph.nodes if node.kind == "record"}

    assert reloaded.version_id == edited.version_id
    assert [event.output_version_id for event in reloaded.history] == [
        event.output_version_id for event in edited.history
    ]
    assert {event.output_version_id for event in edited.history if event.output_version_id} <= record_node_ids
    assert "project_versions:v0" in record_node_ids


def test_project_db_migrates_records_without_version_columns(tmp_path) -> None:
    path = tmp_path / "legacy.sqlite"
    connection = sqlite3.connect(path)
    connection.executescript(
        """
        CREATE TABLE sequences (
            id INTEGER PRIMARY KEY,
            record_id TEXT NOT NULL UNIQUE,
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
        CREATE TABLE features (
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
        CREATE TABLE primers (
            id INTEGER PRIMARY KEY,
            sequence_id INTEGER NOT NULL REFERENCES sequences(id) ON DELETE CASCADE,
            name TEXT NOT NULL,
            sequence TEXT NOT NULL,
            start INTEGER,
            end INTEGER,
            strand INTEGER NOT NULL,
            target_id TEXT
        );
        CREATE TABLE provenance_events (
            id INTEGER PRIMARY KEY,
            sequence_id INTEGER NOT NULL REFERENCES sequences(id) ON DELETE CASCADE,
            event_order INTEGER NOT NULL,
            operation TEXT NOT NULL,
            timestamp TEXT,
            input_record_id TEXT,
            input_record_ids_json TEXT NOT NULL,
            output_record_id TEXT,
            parameters_json TEXT NOT NULL,
            affected_ranges_json TEXT NOT NULL,
            description TEXT NOT NULL
        );
        CREATE TABLE tags (id INTEGER PRIMARY KEY, name TEXT NOT NULL UNIQUE);
        CREATE TABLE sequence_tags (
            sequence_id INTEGER NOT NULL REFERENCES sequences(id) ON DELETE CASCADE,
            tag_id INTEGER NOT NULL REFERENCES tags(id) ON DELETE CASCADE,
            PRIMARY KEY (sequence_id, tag_id)
        );
        CREATE TABLE notes (
            id INTEGER PRIMARY KEY,
            sequence_id INTEGER NOT NULL REFERENCES sequences(id) ON DELETE CASCADE,
            text TEXT NOT NULL,
            created_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
        );
        """
    )
    connection.execute(
        """
        INSERT INTO sequences
            (record_id, name, description, sequence, molecule_type, topology, length, source_path)
        VALUES ('legacy', NULL, NULL, 'ATGC', 'DNA', 'linear', 4, NULL)
        """
    )
    connection.commit()
    connection.close()

    db = ProjectDatabase(path)
    loaded = db.get_record("legacy")

    assert loaded.version_id == "legacy:v0"
    db.update_record("legacy", loaded.insert(2, "GG"))
    assert db.get_record("legacy").history[-1].output_version_id is not None
