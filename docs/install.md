# Installing and Building PlasmidLab

PlasmidLab currently supports local developer installs, CLI usage, and local desktop
bundles with PyInstaller.

## Requirements

- Python 3.12 or newer. Python 3.11 and older are unsupported.
- A working C/C++ runtime suitable for PySide6 and Biopython wheels
- PowerShell for the provided Windows helper scripts

## Editable Install

Use a clean Python 3.12+ virtual environment for local development and distribution
builds. Old global backports of standard-library modules can confuse PyInstaller, so
avoid building from a reused site-packages directory. If your default `python` points
to Python 3.11 or older, invoke a Python 3.12+ interpreter explicitly before creating
the environment.

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python scripts/check_python.py
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

The release scripts also accept `PLASMIDLAB_PYTHON` when your default `python` is not
Python 3.12+:

```powershell
$env:PLASMIDLAB_PYTHON = "C:\Path\To\Python312\python.exe"
```

Launch the desktop GUI:

```powershell
plasmidlab-gui
```

Run the CLI:

```powershell
plasmidlab
plasmidlab digest input.gb --enzymes EcoRI,BamHI
plasmidlab pcr template.gb --fwd ATG... --rev TTA...
plasmidlab convert input.gb output.fasta
plasmidlab map-export input.gb map.svg
```

The CLI commands do not import or require GUI modules unless you run the separate
`plasmidlab-gui` entry point.

## Local Checks

The helper scripts create or reuse an isolated `.venv-build` environment and install
dependencies with `constraints.txt`. This keeps global Python packages out of lint,
test, package, and PyInstaller builds.

```powershell
.\scripts\lint.ps1
.\scripts\test.ps1
```

## Python Package Build

Build source and wheel artifacts under `dist/`:

```powershell
.\scripts\build_package.ps1
```

The script creates a clean Python 3.12+ `.venv-build`, installs PlasmidLab with
development dependencies constrained by `constraints.txt`, and runs:

```powershell
python -m build --no-isolation
```

Bare `python -m build` is also guarded by PlasmidLab's PEP 517 backend wrapper and fails
under Python 3.11 or older. This prevents accidental release metadata builds from an
unsupported interpreter when the helper scripts are bypassed.

## Desktop Bundle With PyInstaller

Generate the placeholder icon and build the desktop app from a clean Python 3.12+
virtual environment:

```powershell
.\scripts\build_pyinstaller.ps1
```

The local desktop artifact is created under:

```text
dist/PlasmidLab/
```

The PyInstaller configuration is in `build_support/pyinstaller/plasmidlab.spec`. It bundles
the PySide6 desktop app, packaged PlasmidLab resources, the built-in open feature
library from `plasmidlab.data.features`, and the original placeholder icon assets.
The spec intentionally uses focused Biopython hidden imports instead of collecting every
`Bio` submodule.

The Python wheel also includes the built-in feature-library JSON files under
`plasmidlab/data/features/`. This keeps Detect Features working from normal wheel
installs as well as editable installs. Lab-specific feature libraries should live outside
the package and be loaded by passing their file or directory path to the feature-library
loader.

## All-In-One Local Build

Run lint, tests, package build, PyInstaller build, post-build smoke checks, and checksum
generation from a clean Python 3.12+ environment:

```powershell
.\scripts\build_all.ps1
```

The release build performs these steps in order:

1. Remove `build/` and `dist/` artifacts.
2. Create a clean `.venv-build` using Python 3.12+.
3. Install constrained development dependencies from `constraints.txt`.
4. Run `pytest`.
5. Run `ruff`.
6. Build wheel and sdist.
7. Install the wheel into an isolated smoke-test venv.
8. Run source-isolated CLI and feature-library smoke checks from the installed wheel.
9. Build the PyInstaller desktop artifact.
10. Run frozen-app smoke checks.
11. Generate SHA256 checksums.

Post-build smoke checks verify:

- wheel install imports `plasmidlab`
- `plasmidlab --help` works from the installed console script
- `plasmidlab digest --help`, `pcr --help`, `convert --help`, and `map-export --help`
  work from the installed console script
- the packaged default feature library loads
- a synthetic FASTA digest workflow runs from the installed wheel
- the PyInstaller executable exists
- the PyInstaller executable can run headless `--help`
- the PyInstaller executable can run headless `--version`
- the PyInstaller executable can run `--smoke-json` and load bundled feature data through
  package resources
- the GUI initializes with Qt offscreen rendering

The installed-wheel smoke checks run from a temporary directory with `PYTHONPATH` and
`PYTHONHOME` removed, so they cannot pass by importing from the source checkout.

The script writes `dist/SHA256SUMS.txt` for the wheel, sdist, and zipped PyInstaller
artifact. Checksum generation fails if the wheel, sdist, or PyInstaller directory is
missing.

## Release Dependency Constraints

`constraints.txt` pins the v0.1-alpha release-build dependency set. Refresh the lock
metadata only from a clean Python 3.12+ environment:

```powershell
.\scripts\lock_requirements.ps1
```

The command writes `requirements.lock`; review it and update `constraints.txt` before
tagging a release.

Do not tag or publish v0.1-alpha until `scripts/build_all.ps1` passes on Python 3.12+
or Python 3.13.

## Continuous Integration

`.github/workflows/ci.yml` runs on Windows with Python 3.12 and 3.13 and invokes
`scripts/build_all.ps1`. CI therefore uses the same clean-venv, constrained dependency,
smoke-test, PyInstaller, and checksum path as local release builds. Artifact uploads are
configured to fail if the expected wheel, sdist, PyInstaller directory, zipped desktop
artifact, or checksum file is missing. The workflow invokes the PowerShell build script
with process execution-policy bypass so Windows runner policy does not block the release
validation script.

When a maintainer does not have local Python 3.12+, the GitHub Actions workflow is the
authoritative release-validation path. v0.1-alpha remains `NO-GO` until the Python 3.12+
CI job passes end to end and uploads the package artifacts, desktop artifact, and
`SHA256SUMS.txt`.
