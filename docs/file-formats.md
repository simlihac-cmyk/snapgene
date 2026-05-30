# File Formats

PlasmidLab uses open interchange formats. It does not read or write proprietary
SnapGene `.dna` files.

## FASTA

FASTA import/export is implemented in `plasmidlab.io.fasta`.

Imported fields:

- Sequence text.
- Record id.
- Record name when available from Biopython.
- Description line.

FASTA does not carry standard feature annotations, molecule topology, or rich project
metadata. PlasmidLab therefore defaults imported FASTA records to DNA and linear
topology unless the caller provides explicit `molecule_type` or `topology` arguments.

## GenBank

GenBank import/export is implemented in `plasmidlab.io.genbank` using Biopython parsing
and writing primitives.

Imported fields:

- Sequence text.
- Record id, name, and description.
- Molecule type from GenBank annotations when available.
- Linear or circular topology from GenBank annotations when available.
- Feature type, zero-based half-open feature segments, strand, and qualifiers.
- Repeated qualifier values as ordered lists, represented internally as immutable tuples.
- Compound location operators for supported `join(...)` and `order(...)` locations.

Exported fields:

- LOCUS metadata with molecule type, topology, synthetic division, and deterministic
  placeholder date.
- Sequence text.
- Features as GenBank simple or compound locations.
- Qualifiers stored on PlasmidLab features.
- Repeated qualifier values as repeated GenBank qualifiers.
- Supported compound location operators (`join` and `order`) when present in imported
  metadata or explicitly set by PlasmidLab.

Internally, all imported GenBank locations are normalized to zero-based half-open
coordinates. On export, Biopython converts those internal coordinates back to GenBank's
one-based inclusive location display.

Compound GenBank locations are stored as multiple `FeatureSegment` values. Circular
features that cross the origin are represented as ordered segments, for example
`[112, 120)` followed by `[0, 6)` on a 120 bp circular plasmid.

PlasmidLab preserves segment order from Biopython for compound locations. This matters
for reverse-strand CDS records, where a valid `complement(join(...))` feature may be
stored in transcript order rather than ascending coordinate order.

### GenBank Loss Warnings

`plasmidlab.io.genbank` emits `GenBankLossyWarning` when a GenBank construct can be read
only with explicit degradation. The warning text is also kept on the imported feature as
PlasmidLab metadata under `_plasmidlab_location_warnings` so downstream code can surface
the issue to users.

Supported without warning:

- Simple local feature locations.
- Compound local `join(...)` locations.
- Compound local `order(...)` locations.
- Forward, reverse, and unknown feature strands.
- Multi-valued qualifiers, including repeated `note`, `db_xref`, and `translation`
  qualifiers.
- `source` features and LOCUS topology metadata.

Imported with explicit warning:

- Fuzzy boundaries such as `<3..>8`. PlasmidLab stores the numeric zero-based half-open
  boundary reported by Biopython and records warning metadata. Export uses exact numeric
  boundaries because the core model does not yet represent fuzzy boundary classes.

Rejected with a clear error:

- Remote feature locations or cross-record references such as `J00194.1:1..10`.
  PlasmidLab does not silently convert these into local coordinates.

Internal PlasmidLab metadata qualifiers beginning with `_plasmidlab_` are used to retain
location operator and warning information. These metadata keys are not written as GenBank
qualifiers; they are translated back into supported GenBank locations where possible.
