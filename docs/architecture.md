# PlasmidLab Architecture

PlasmidLab separates biological logic, file IO, rendering, GUI concerns, and project
history so that core behavior can be tested deterministically without a desktop runtime.

## Package Boundaries

- `plasmidlab.core` contains pure sequence records, annotations, primers, enzyme sites,
  transformations, and public molecular biology algorithms.
- `plasmidlab.io` owns open file-format import/export such as FASTA, GenBank, EMBL, and
  PlasmidLab JSON project files.
- `plasmidlab.render` converts records and annotations into display-independent drawing
  primitives and export formats. Map overlays from restriction analyses and primer
  binding searches are prepared here so GUI code does not fabricate biological
  annotations for display.
- `plasmidlab.gui` will contain PySide6 widgets, windows, actions, and application state.
  GUI code may call core services but must not implement biological algorithms.
- `plasmidlab.history` records provenance events and undo/redo stacks.

## Coordinate Conventions

All sequence coordinates are stored as zero-based, half-open intervals: `[start, end)`.
The first base is index `0`; an interval ending at `n` excludes base `n`. This matches
normal Python slicing and avoids ambiguity when editing sequence ranges.

Examples:

- `[0, 1)` selects the first base.
- `[5, 10)` selects five bases: indices `5`, `6`, `7`, `8`, and `9`.
- `[3, 3)` is an empty interval and is not valid for biological features.

External file formats may use other conventions. Importers must normalize coordinates to
zero-based half-open intervals on read, and exporters must convert back to the target
format on write.

The GUI may convert coordinates into screen positions for highlighting or table display,
but biological coordinate construction, such as map overlays for enzyme sites and primer
bindings, belongs in core or render APIs.

## Circular Sequence Handling

Records carry an explicit topology: `linear` or `circular`. Algorithms must not infer
circularity from sequence length or feature shape.

For storage, feature intervals are still represented in normalized zero-based half-open
coordinates where `0 <= start < end <= sequence_length`. Features that cross the origin
of a circular molecule should be represented by multiple intervals or by a future
compound-location model rather than by negative indices or `end > length`.

Circular algorithms, such as primer binding across the origin or restriction digest on a
plasmid, may use doubled-sequence search internally. Returned coordinates must be
normalized back to the record coordinate system.

Circular interval APIs must distinguish zero-length selections from full-circle
selections. In general sequence-editing APIs, `start == end` means a zero-length
half-open interval. A full-circle interval must be represented explicitly by the
operation result or parameters rather than inferred from equal coordinates.

`SequenceRecord.replace(start, end, sequence)` preserves the stored origin. Circular
replacements where `start <= end` behave like normal half-open replacement without
rotating the sequence. Origin-crossing circular replacement where `start > end` is
rejected because replacing two ranges on either side of the origin would otherwise
require either moving the origin or inventing an implicit cut. Callers should linearize
at an explicit cut or use delete/insert steps with clear coordinates.

Circular PCR is one of the APIs that can explicitly report a full-circle interval. A PCR
product with equal start and end coordinates may be marked `full_circle=True`; its
provenance affected ranges are reported as `[(start, length), (0, start)]`, or
`[(0, length)]` when `start == 0`, never as a zero-length interval.

PCR with tailed primers treats the 3-prime suffix as the annealing region and the
5-prime prefix as a non-template tail. Binding search considers a configurable
annealing-length range instead of stopping at the first hit. Candidate bindings are
ranked by annealing length, mismatch count, Tm, uniqueness, and optional target-region
fit; multiple plausible products remain explicit ambiguity unless the caller selects a
binding pair or provides a target region that identifies one product.

Restriction digest, restriction cloning, and Golden Gate assembly use explicit physical
fragment ends derived from top-strand and bottom-strand cut geometry. See
`docs/restriction-overhangs.md` for the overhang sequence, left/right end, Type IIS, and
ambiguous-base compatibility rules.

## Transformations and Provenance

Sequence transformations should be pure from the caller's perspective: return a new
sequence record and a provenance event instead of mutating a record in place. Events
should capture the operation name, inputs, parameters, timestamp, and enough context to
support project history and undo/redo.

Annotation preservation is part of the transformation contract. Insertions, deletions,
reverse complements, and circular slicing should keep features when coordinates can be
mapped unambiguously and should record when annotations are truncated or dropped.

Feature annotations are user-authored biological intervals, so edit operations may
remap, truncate, or drop them using deterministic coordinate projection rules. Each edit
stores an `annotation_semantics` provenance payload with feature counts before and after
the operation.

Primer entries are treated as primer definitions plus optional binding coordinates. Edits
that alter sequence coordinates invalidate existing primer binding coordinates rather
than shifting them blindly; the primer name and sequence are retained with an unknown
binding location. Reverse complement operations may retain primer binding annotations by
mapping `[start, end)` to `[length - end, length - start)` and flipping strand.

Enzyme-site annotations are derived from a selected enzyme set. Sequence edits and
topology-changing operations drop stored enzyme-site coordinates and record that the
sites must be recomputed from the edited sequence and the selected enzyme set. They must
not be shifted silently.

## Legal Boundaries

PlasmidLab is an independent application. It must not copy SnapGene code, UI layouts,
icons, visual assets, proprietary datasets, reverse-engineered proprietary file formats,
or distinctive product presentation.

The project may implement common molecular biology workflows using public literature,
open specifications, and independently written code. Supported interchange formats should
be open formats such as FASTA, GenBank, EMBL, SVG, PDF, CSV/TSV, and PlasmidLab's own
documented JSON project format.
