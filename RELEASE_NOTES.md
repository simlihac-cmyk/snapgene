# Release Notes

## v0.1.0-alpha

### Summary

PlasmidLab v0.1.0-alpha is the first alpha release candidate for local molecular
biology research and evaluation workflows. It is an independent open-format desktop
application project built from original code, public molecular biology algorithms, and
open data definitions.

This alpha is not a drop-in replacement for any commercial molecular biology package.
It is intended for local evaluation, reproducibility testing, and internal research
workflows.

### Implemented Capabilities

- Immutable/copy-on-write core sequence records, feature segments, primers, enzyme
  sites, provenance events, molecule topology, and molecule type models.
- Zero-based, half-open coordinate conventions with explicit circular and linear
  molecule handling.
- FASTA and GenBank import/export, including common GenBank feature, qualifier,
  topology, compound-location, and warning behavior.
- Restriction enzyme search, digest simulation, sequence-derived physical overhang
  modeling, and restriction cloning support.
- Primer binding, primer metrics, primer3-backed primer design, and PCR simulation with
  circular templates and primer tails.
- ORF detection, translation, CDS validation, and reverse-translation utilities.
- Gibson-like overlap assembly, Golden Gate Type IIS assembly, inverse PCR, and
  mutagenesis simulations using synthetic/open examples.
- Plasmid map drawing model, SVG export, feature/enzyme/primer overlays, and a PySide6
  map widget.
- PySide6 desktop GUI with map, sequence, features, enzymes, primers, gel, history,
  project search, and Sanger trace views.
- Pairwise alignment, reference discrepancy reports, and AB1 Sanger trace import/viewing
  based on called bases.
- Approximate agarose gel simulation with open JSON ladder definitions and SVG/PNG
  export paths.
- Open JSON feature-library loading through package resources plus optional user
  libraries.
- History/provenance tracking, undo/redo support, history graph export, SQLite project
  collections, search, batch operations, CLI entry points, and PyInstaller packaging.

### Validation Results

Release validation was performed by the GitHub Actions workflow `CI`.

- Python 3.12 release-validation: passed.
- Python 3.13 release-validation: passed.
- `pytest`: passed, expected at 150 tests or the CI-reported value.
- `ruff check src tests scripts`: passed.
- `python -m build`: passed.
- Installed-wheel smoke: passed.
- Feature-library loading through `importlib.resources`: passed.
- CLI `plasmidlab --help`: passed.
- CLI subcommand help for `digest`, `pcr`, `convert`, and `map-export`: passed.
- Synthetic FASTA digest smoke workflow: passed.
- PyInstaller build: passed.
- Frozen executable smoke: `--help`, `--version`, and `--smoke-json` passed.
- SHA256 checksum generation: passed.

### Artifacts

The passing `CI` workflow uploaded these artifact groups:

- `PlasmidLab-windows-py3.12`
- `python-package-py3.12`
- `PlasmidLab-windows-py3.13`
- `python-package-py3.13`

The Python package artifacts include the wheel, sdist, and `SHA256SUMS.txt`. The desktop
artifacts include the PyInstaller directory, zipped Windows desktop artifact, and
`SHA256SUMS.txt`.

### Python Requirement

PlasmidLab v0.1.0-alpha requires Python 3.12 or newer. Python 3.11 and older are
unsupported. The build backend and release scripts fail early when run with unsupported
Python versions.

### Known Limitations

- PlasmidLab v0.1.0-alpha is intended for local research/evaluation workflows.
- It is not validated for clinical, diagnostic, regulated, or production
  decision-making use.
- GenBank import/export supports common constructs; unsupported or lossy constructs may
  warn or fail explicitly.
- Golden Gate and restriction cloning use sequence-derived physical overhangs, but
  critical designs should still be experimentally verified.
- Sanger AB1 support uses base calls present in trace files and does not perform de novo
  base calling.
- Agarose gel simulation is approximate visualization, not an empirical electrophoresis
  predictor.
- Proprietary SnapGene `.dna` import/export is not supported.

### Unsupported Proprietary Formats

PlasmidLab does not support proprietary SnapGene `.dna` import/export and does not copy
SnapGene code, UI assets, proprietary data, feature databases, or proprietary file-format
details.

### Scientific Assumptions

- Core coordinates are zero-based, half-open intervals.
- Circular and linear molecules are represented explicitly.
- Restriction and Golden Gate overhang compatibility is based on physical fragment-end
  geometry and sequence-derived overhangs.
- Primer/PCR simulations model deterministic sequence-level behavior and should not be
  treated as experimental validation.
- Agarose gel output is an approximate visualization model.

### Checksum References

Use the `SHA256SUMS.txt` file uploaded with each artifact group to verify the wheel,
sdist, and zipped PyInstaller desktop artifact before local evaluation.
