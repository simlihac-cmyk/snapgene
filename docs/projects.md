# Project Collections and Search

PlasmidLab provides a SQLite-backed project collection layer in
`plasmidlab.project.ProjectDatabase`.

## Database Contents

The project database stores normalized biological and project metadata:

- `sequences`: record ID, name, description, sequence, molecule type, topology, length
- `features`: feature name/type, strand, zero-based ranges, qualifiers
- `primers`: primer name, sequence, optional binding interval
- `provenance_events`: enriched event history for each record
- `tags`: user tags with a sequence/tag join table
- `notes`: user notes

Searchable text is indexed from sequence names/descriptions, feature names, primer
names, and notes. SQLite FTS5 is used when available, with a LIKE fallback for SQLite
builds without FTS5.

## Text Search

Use `ProjectDatabase.search_text(query)` to search:

- sequence IDs, names, and descriptions
- feature names
- primer names
- notes

The result reports the record primary key, record ID, matched field, and matched text.

## Sequence Search

Use `ProjectDatabase.search_sequence(query, mismatches=0)` for exact or approximate
sequence search. DNA records are searched on the forward strand and, by default, the
reverse-complement strand. Circular records support origin-crossing matches.

Protein records are searched as ordinary exact or mismatch-tolerant strings without
reverse-complement handling.

## Batch Operations

The project database includes batch helpers:

- `import_folder(path)`: import GenBank and FASTA files from a folder
- `export_records(records, output_dir, file_format=...)`: export selected records
- `batch_detect_features(records)`: run automatic feature annotation
- `batch_restriction_analysis(records, enzymes)`: run restriction analysis
- `batch_export_map_svg(records, output_dir)`: export selected maps as SVG

The GUI project panel can open a SQLite project database, import a folder into a
project, and search/filter visible records.
