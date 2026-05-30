# AGENTS.md

Guidance for coding agents working on PlasmidLab.

## Project Mission

Build an independent molecular biology desktop app for plasmid visualization,
annotation, restriction analysis, primer design, PCR and cloning simulation, sequence
editing, ORF/translation, alignment, Sanger trace verification, gel simulation, and
project history.

Do not copy SnapGene code, UI, assets, proprietary data, or proprietary file formats.
Use original implementation work, public algorithms, open specifications, and synthetic
example data.

## Repository Layout

- `src/plasmidlab/core/`: deterministic biological data models and algorithms.
- `src/plasmidlab/io/`: import/export for FASTA, GenBank, EMBL, and JSON projects.
- `src/plasmidlab/gui/`: PySide6 GUI code only.
- `src/plasmidlab/render/`: plasmid map and sequence rendering.
- `src/plasmidlab/history/`: provenance and undo/redo support.
- `tests/`: unit and integration tests.
- `examples/`: synthetic example plasmids only.
- `docs/`: architecture and contributor-facing design notes.

## Conventions

- Store all sequence coordinates as zero-based, half-open intervals: `[start, end)`.
- Model circular and linear molecules explicitly.
- Keep core algorithms deterministic and free of GUI imports.
- Do not implement biological logic in GUI modules.
- Each sequence transformation should return a new record plus a provenance event.
- Preserve annotations across editing operations whenever possible.
- Add tests with or before implementation changes.

## Development Commands

```bash
python -m pip install -e ".[dev]"
python -m pytest
```

Prefer focused tests for new core behavior. Broaden integration tests when changes cross
IO, history, rendering, or GUI boundaries.
