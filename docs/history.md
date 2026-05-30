# History, Provenance, Undo, and Redo

PlasmidLab records biological transformations as `ProvenanceEvent` objects on
`SequenceRecord.history`. Events are stored with records and remain pure data, so they
can be tested without the GUI.

## Event Schema

Each new event includes:

- `operation`: stable operation name, such as `insert`, `delete`, `pcr`, or
  `gibson_assembly`
- `timestamp`: UTC ISO timestamp
- `input_record_ids`: one or more source record IDs
- `output_record_id`: derived record or product ID
- `input_version_ids`: one or more source sequence-version IDs
- `output_version_id`: derived sequence-version ID
- `parameters`: operation-specific structured parameters
- `affected_ranges`: zero-based half-open coordinate ranges affected by the operation
- `description`: short human-readable summary

The older `input_record_id` field remains available as a compatibility alias for the
first input record.

`SequenceRecord.id` is the stable biological record or molecule identifier. Each
`SequenceRecord` also has a `version_id`, which identifies one concrete sequence version.
Core edits keep the stable `id` and create a new `version_id`; product-producing
workflows such as PCR and cloning create product record IDs with their own output version
IDs. Older projects without version fields load with migration-safe fallback version IDs.

## Core Coverage

Core sequence edits add provenance directly:

- reverse complement
- circularize and linearize
- slice
- insert
- delete
- replace

Cloning and amplification workflows also create provenance:

- restriction cloning
- Gibson-like assembly
- Golden Gate assembly
- inverse-PCR mutagenesis
- PCR products

PCR still returns `PCRProduct`, but the product now carries a `provenance` event.

## History Graph

`plasmidlab.history.build_history_graph()` converts records and events into a graph:

- record nodes represent concrete sequence versions, keyed by `version_id`
- event nodes represent operations
- input edges connect source versions to events
- output edges connect events to derived output versions

Use `write_history_json(records, path)` to export the graph as JSON.

## GUI Undo/Redo

The desktop GUI stores immutable `SequenceRecord` snapshots for editing undo/redo.
Undo restores the exact previous sequence and annotations. Redo restores the exact
edited sequence and annotations. The History tab shows both an event table and a simple
graph visualization, and the File menu can export history JSON.
