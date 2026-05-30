# PlasmidLab

PlasmidLab is an independent molecular biology desktop app project for plasmid and
sequence workflows. The long-term goal is a SnapGene-like research workflow tool built
from original code, public molecular biology algorithms, and open file formats.

This repository is an early implementation with deterministic core algorithms, open
file-format IO, a GUI-independent plasmid-map renderer, and a first usable PySide6
desktop shell.

## Implemented

- Python `src/` package layout.
- Initial package metadata in `pyproject.toml`.
- Core package namespace with immutable sequence, feature, primer, enzyme-site, and
  provenance models.
- Zero-based, half-open interval validation for core features.
- Compound feature segments, strand support, explicit topology, and molecule-type enums.
- Core sequence edits: reverse complement, circularize, linearize, slice, insert,
  delete, and replace.
- FASTA import/export.
- GenBank import/export with feature segments, strands, topology, and qualifiers.
- Restriction enzyme site search and digest simulation for linear and circular DNA,
  using Biopython enzyme metadata.
- Primer binding search, primer metrics, primer3-backed primer design wrapper, and
  PCR simulation with circular templates and 5-prime primer tails.
- ORF detection, translation, CDS validation, and reverse translation utilities.
- Core cloning simulations for restriction cloning, Gibson-like overlap assembly,
  Golden Gate Type IIS assembly, and inverse-PCR mutagenesis.
- Pairwise global, local, and semi-global alignment for DNA/protein-style sequences,
  reference discrepancy reports, and TSV report export.
- Sanger AB1 import using Biopython, called-base/quality/chromatogram extraction,
  trace-to-reference alignment, and a PySide6 chromatogram widget.
- Approximate agarose gel simulation for digest/PCR fragments, open JSON ladder
  definitions, SVG export, and PySide6 PNG rendering.
- Automatic common feature detection from packaged open JSON libraries under
  `plasmidlab.data.features`, plus user-provided JSON library paths, with
  exact/near-exact matching and GUI review before applying.
- Enriched provenance events, history graph JSON export, and GUI undo/redo for record
  edits.
- SQLite-backed project collections with sequences, features, primers, provenance,
  tags, notes, full-text search, sequence search, and batch import/export/analysis.
- GUI-independent plasmid map drawing model with SVG export, plus a minimal PySide6
  widget that renders the drawing model.
- First usable PySide6 desktop shell with open/save/export, map, sequence, features,
  enzymes, primers, gel, history tabs, project search panel, and standalone AB1 trace
  windows.
- CLI entry points for digest, PCR, format conversion, and SVG map export.
- Local distribution tooling with PyInstaller configuration, original placeholder icon
  assets, constrained Python 3.12+ PowerShell build scripts, release smoke checks,
  checksums, and GitHub Actions CI.

Launch the desktop GUI:

```bash
plasmidlab-gui
```
- CLI digest command:

```bash
plasmidlab digest input.gb --enzymes EcoRI,BamHI
```

- CLI PCR command:

```bash
plasmidlab pcr template.gb --fwd ATG... --rev TTA...
```

- CLI convert and map export commands:

```bash
plasmidlab convert input.gb output.fasta
plasmidlab map-export input.gb map.svg
```

- Placeholder packages for IO, rendering, GUI, history, examples, and tests.
- Minimal pytest suite proving package import and coordinate validation.
- Architecture notes in `docs/architecture.md`.
- Local install and desktop build notes in `docs/install.md`.

## Planned

- EMBL and JSON project import/export.
- Rich alignment review and Sanger discrepancy navigation in the main project window.
- Exportable reports around maps, gels, alignments, and cloning workflows.
- PDF map export, signed installers, and richer release packaging.

## Development

PlasmidLab requires Python 3.12 or newer. Python 3.11 and older are unsupported.
Create a clean Python 3.12+ environment, install the package with development extras,
and run the tests:

```bash
python scripts/check_python.py
python -m pip install -e ".[dev]"
python -m pytest
```

Release builds should use the provided clean-venv scripts so global Python packages do
not affect packaging:

```bash
pwsh ./scripts/build_all.ps1
```

`constraints.txt` defines the v0.1-alpha release-build dependency set. Refresh
`requirements.lock` and review the constraints only from a clean Python 3.12+
environment.

The core package must remain deterministic and testable without importing GUI code.
GUI modules will use PySide6 later, but biological logic belongs under
`src/plasmidlab/core/` and file-format logic belongs under `src/plasmidlab/io/`.

## Legal Boundaries

PlasmidLab is not affiliated with SnapGene. The project must not copy SnapGene code,
user interface assets, proprietary data, reverse-engineered proprietary formats, or
distinctive product presentation. Features should be implemented independently from
public scientific methods and open specifications. Feature libraries in this repository
are open synthetic examples, not copied commercial databases.
