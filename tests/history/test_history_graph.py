import json

from plasmidlab.core import Feature, SequenceRecord
from plasmidlab.history import build_history_graph, history_graph_to_json, write_history_json


def test_sequence_edits_record_enriched_provenance_and_graph_export(tmp_path) -> None:
    record = SequenceRecord(
        id="history",
        sequence="AAAA",
        features=(Feature(type="misc_feature", start=1, end=3, name="feature"),),
    )
    edited = record.insert(2, "GG").delete(0, 1)

    insert_event = edited.history[0]
    delete_event = edited.history[1]
    graph = build_history_graph(edited)
    path = tmp_path / "history.json"
    write_history_json(edited, path)
    exported = json.loads(path.read_text(encoding="utf-8"))

    assert insert_event.operation == "insert"
    assert insert_event.timestamp
    assert insert_event.input_record_ids == ("history",)
    assert insert_event.input_version_ids == ("history:v0",)
    assert insert_event.output_record_id == "history"
    assert insert_event.output_version_id != "history:v0"
    assert delete_event.input_version_ids == (insert_event.output_version_id,)
    assert delete_event.output_version_id == edited.version_id
    assert insert_event.affected_ranges == ((2, 4),)
    assert "Inserted 2 bases" in insert_event.description
    assert delete_event.affected_ranges == ((0, 1),)
    assert any(node.kind == "event" and node.label == "insert" for node in graph.nodes)
    assert any(edge.label == "output" for edge in graph.edges)
    assert exported["nodes"]
    assert exported["edges"]
    assert history_graph_to_json(graph).startswith("{")


def test_history_graph_uses_distinct_version_nodes_for_sequential_edits() -> None:
    original = SequenceRecord(id="versions", sequence="AAAA")
    inserted = original.insert(1, "GG")
    deleted = inserted.delete(0, 1)
    graph = build_history_graph(deleted)

    record_nodes = {node.id for node in graph.nodes if node.kind == "record"}

    assert original.version_id == "versions:v0"
    assert inserted.version_id != original.version_id
    assert deleted.version_id != inserted.version_id
    assert {original.version_id, inserted.version_id, deleted.version_id} <= record_nodes
    assert len(record_nodes) == 3
    assert any(edge.source == original.version_id and edge.label == "input" for edge in graph.edges)
    assert any(edge.target == deleted.version_id and edge.label == "output" for edge in graph.edges)
