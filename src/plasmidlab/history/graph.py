"""History graph construction and JSON export."""

from __future__ import annotations

import json
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any

from plasmidlab.core.models import ProvenanceEvent, SequenceRecord


@dataclass(frozen=True, slots=True)
class HistoryNode:
    """A sequence-record or operation-event node."""

    id: str
    kind: str
    label: str
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))


@dataclass(frozen=True, slots=True)
class HistoryEdge:
    """A derivation edge between records and operation events."""

    source: str
    target: str
    label: str = ""


@dataclass(frozen=True, slots=True)
class HistoryGraph:
    """A lightweight provenance graph."""

    nodes: tuple[HistoryNode, ...]
    edges: tuple[HistoryEdge, ...]


def build_history_graph(records: SequenceRecord | Iterable[SequenceRecord]) -> HistoryGraph:
    """Build a graph where record nodes connect through operation-event nodes."""

    record_tuple = (records,) if isinstance(records, SequenceRecord) else tuple(records)
    nodes: dict[str, HistoryNode] = {}
    edges: list[HistoryEdge] = []
    for record in record_tuple:
        _ensure_version_node(
            nodes,
            record.version_id,
            record_id=record.id,
            sequence_length=record.length if not record.history else None,
        )
        previous_versions: dict[str, str] = {record.id: _initial_version_id(record.id)}
        for index, event in enumerate(record.history):
            output_record_id = event.output_record_id or record.id
            output_version_id = event.output_version_id or _legacy_output_version_id(
                output_record_id,
                index,
                record,
            )
            event_id = _event_node_id(output_version_id, index, event)
            nodes[event_id] = HistoryNode(
                id=event_id,
                kind="event",
                label=event.operation,
                metadata={
                    "timestamp": event.timestamp,
                    "operation": event.operation,
                    "parameters": dict(event.parameters),
                    "affected_ranges": event.affected_ranges,
                    "description": event.description,
                    "input_record_ids": event.input_record_ids,
                    "input_version_ids": event.input_version_ids,
                    "output_record_id": output_record_id,
                    "output_version_id": output_version_id,
                },
            )
            input_record_ids = event.input_record_ids or (
                (event.input_record_id,) if event.input_record_id else ()
            )
            input_version_ids = event.input_version_ids or (
                (event.input_version_id,) if event.input_version_id else ()
            )
            if input_version_ids:
                inputs = tuple(
                    (
                        input_version_ids[position],
                        input_record_ids[position] if position < len(input_record_ids) else None,
                    )
                    for position in range(len(input_version_ids))
                )
            else:
                inputs = tuple(
                    (
                        previous_versions.get(input_id, _initial_version_id(input_id)),
                        input_id,
                    )
                    for input_id in input_record_ids
                )
            for input_version_id, input_record_id in inputs:
                _ensure_version_node(nodes, input_version_id, record_id=input_record_id)
                edges.append(HistoryEdge(source=input_version_id, target=event_id, label="input"))
            previous_versions[output_record_id] = output_version_id
            _ensure_version_node(
                nodes,
                output_version_id,
                record_id=output_record_id,
                sequence_length=record.length if index == len(record.history) - 1 else None,
            )
            edges.append(HistoryEdge(source=event_id, target=output_version_id, label="output"))
    return HistoryGraph(nodes=tuple(nodes.values()), edges=tuple(_dedupe_edges(edges)))


def history_graph_to_dict(graph: HistoryGraph) -> dict[str, Any]:
    """Convert a history graph into JSON-safe data."""

    return {
        "nodes": [
            {
                "id": node.id,
                "kind": node.kind,
                "label": node.label,
                "metadata": _json_safe(dict(node.metadata)),
            }
            for node in graph.nodes
        ],
        "edges": [
            {"source": edge.source, "target": edge.target, "label": edge.label}
            for edge in graph.edges
        ],
    }


def history_graph_to_json(graph: HistoryGraph, *, indent: int = 2) -> str:
    """Serialize a history graph as JSON text."""

    return json.dumps(history_graph_to_dict(graph), indent=indent, sort_keys=True)


def write_history_json(
    records: SequenceRecord | Iterable[SequenceRecord],
    path: str | Path,
    *,
    indent: int = 2,
) -> None:
    """Export record history as JSON."""

    graph = build_history_graph(records)
    Path(path).write_text(history_graph_to_json(graph, indent=indent), encoding="utf-8")


def _ensure_version_node(
    nodes: dict[str, HistoryNode],
    version_id: str,
    *,
    record_id: str | None = None,
    sequence_length: int | None = None,
) -> None:
    if version_id in nodes:
        if sequence_length is not None and "sequence_length" not in nodes[version_id].metadata:
            metadata = dict(nodes[version_id].metadata)
            metadata["sequence_length"] = sequence_length
            nodes[version_id] = HistoryNode(
                id=version_id,
                kind="record",
                label=nodes[version_id].label,
                metadata=metadata,
            )
        return
    metadata: dict[str, Any] = {"version_id": version_id}
    if record_id is not None:
        metadata["record_id"] = record_id
    if sequence_length is not None:
        metadata["sequence_length"] = sequence_length
    nodes[version_id] = HistoryNode(
        id=version_id,
        kind="record",
        label=record_id or version_id,
        metadata=metadata,
    )


def _event_node_id(output_version_id: str, index: int, event: ProvenanceEvent) -> str:
    return f"event:{output_version_id}:{index}:{event.operation}"


def _initial_version_id(record_id: str) -> str:
    return f"{record_id}:v0"


def _legacy_output_version_id(record_id: str, index: int, record: SequenceRecord) -> str:
    if index == len(record.history) - 1:
        return record.version_id
    return f"{record_id}:legacy-v{index + 1}"


def _dedupe_edges(edges: Iterable[HistoryEdge]) -> tuple[HistoryEdge, ...]:
    seen: set[tuple[str, str, str]] = set()
    deduped: list[HistoryEdge] = []
    for edge in edges:
        key = (edge.source, edge.target, edge.label)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(edge)
    return tuple(deduped)


def _json_safe(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, tuple | list):
        return [_json_safe(item) for item in value]
    return value
