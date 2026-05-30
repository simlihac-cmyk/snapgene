"""Post-build smoke checks for PlasmidLab release artifacts."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

MIN_VERSION = (3, 12)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dist", type=Path, default=Path("dist"))
    parser.add_argument("--constraints", type=Path, default=Path("constraints.txt"))
    parser.add_argument("--pyinstaller-dir", type=Path, default=Path("dist/PlasmidLab"))
    parser.add_argument("--wheel-only", action="store_true", help="only run wheel install smoke checks")
    parser.add_argument(
        "--pyinstaller-only",
        action="store_true",
        help="only run frozen PyInstaller smoke checks",
    )
    parser.add_argument(
        "--skip-gui-init",
        action="store_true",
        help="skip optional offscreen GUI initialization in the wheel smoke",
    )
    args = parser.parse_args()
    try:
        _check_python_version()
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    if args.wheel_only and args.pyinstaller_only:
        raise ValueError("--wheel-only and --pyinstaller-only are mutually exclusive")

    if not args.pyinstaller_only:
        _smoke_wheel(args.dist, args.constraints, skip_gui_init=args.skip_gui_init)
    if not args.wheel_only:
        _smoke_pyinstaller(args.pyinstaller_dir)

    print("release smoke checks passed")
    return 0


def _check_python_version() -> None:
    if sys.version_info[:2] >= MIN_VERSION:
        return
    detected = ".".join(str(part) for part in sys.version_info[:3])
    msg = f"release smoke checks require Python 3.12 or newer; detected {detected}"
    raise RuntimeError(msg)


def _smoke_wheel(dist: Path, constraints: Path, *, skip_gui_init: bool) -> None:
    wheel = _single_latest(dist.glob("*.whl"), "wheel").resolve()
    _single_latest(dist.glob("*.tar.gz"), "sdist")
    constraints = constraints.resolve()
    isolated_env = _isolated_python_env()
    with tempfile.TemporaryDirectory(prefix="plasmidlab-wheel-smoke-") as temp_dir:
        temp_path = Path(temp_dir)
        venv = temp_path / "venv"
        _run((sys.executable, "-m", "venv", str(venv)), env=isolated_env)
        python = _venv_python(venv)
        _run((str(python), "-m", "pip", "install", "--upgrade", "pip"), env=isolated_env)
        _run(
            (
                str(python),
                "-m",
                "pip",
                "install",
                "--constraint",
                str(constraints),
                str(wheel),
            ),
            env=isolated_env,
        )
        _run((str(python), "-c", _IMPORT_AND_FEATURE_LIBRARY_CHECK), env=isolated_env, cwd=temp_path)
        cli = _console_script(venv, "plasmidlab")
        _run((cli, "--help"), env=isolated_env, cwd=temp_path)
        for subcommand in ("digest", "pcr", "convert", "map-export"):
            _run((cli, subcommand, "--help"), env=isolated_env, cwd=temp_path)
        _run_synthetic_cli_workflow(cli, temp_path, env=isolated_env)
        if not skip_gui_init:
            _run(
                (str(python), "-c", _GUI_INIT_CHECK),
                env=_isolated_python_env({"QT_QPA_PLATFORM": "offscreen"}),
                cwd=temp_path,
            )


def _smoke_pyinstaller(pyinstaller_dir: Path) -> None:
    executable = _pyinstaller_executable(pyinstaller_dir).resolve()
    isolated_env = _isolated_python_env()
    with tempfile.TemporaryDirectory(prefix="plasmidlab-frozen-smoke-") as temp_dir:
        temp_path = Path(temp_dir)
        _run((str(executable), "--help"), env=isolated_env, cwd=temp_path, timeout=30)
        _run((str(executable), "--version"), env=isolated_env, cwd=temp_path, timeout=30)
        smoke_path = temp_path / "frozen-smoke.json"
        _run(
            (str(executable), "--smoke-json", str(smoke_path)),
            env=isolated_env,
            cwd=temp_path,
            timeout=30,
        )
        data = json.loads(smoke_path.read_text(encoding="utf-8"))
        if not data.get("feature_library_entries"):
            raise RuntimeError("frozen smoke did not load bundled feature-library entries")
        if not data.get("feature_library"):
            raise RuntimeError("frozen smoke did not report feature-library availability")
        if not data.get("package_resources"):
            raise RuntimeError("frozen smoke did not report package-resource availability")
        if not data.get("version"):
            raise RuntimeError("frozen smoke did not report PlasmidLab version")


_IMPORT_AND_FEATURE_LIBRARY_CHECK = """
import plasmidlab
from plasmidlab.core.feature_annotation import load_feature_library
library = load_feature_library()
assert plasmidlab.__version__
assert library
print(f"imported plasmidlab {plasmidlab.__version__}; loaded {len(library)} feature entries")
"""


_GUI_INIT_CHECK = """
from PySide6.QtWidgets import QApplication
from plasmidlab.gui.main_window import PlasmidLabMainWindow
app = QApplication.instance() or QApplication([])
window = PlasmidLabMainWindow()
assert window.windowTitle() == "PlasmidLab"
window.close()
print("GUI initialized")
"""


def _single_latest(paths: object, label: str) -> Path:
    path_list = sorted(paths, key=lambda path: path.stat().st_mtime)
    if not path_list:
        raise FileNotFoundError(f"No {label} artifact found")
    return path_list[-1]


def _pyinstaller_executable(path: Path) -> Path:
    executable = path / ("PlasmidLab.exe" if os.name == "nt" else "PlasmidLab")
    if not executable.exists():
        raise FileNotFoundError(f"PyInstaller executable not found: {executable}")
    print(f"found PyInstaller executable: {executable}")
    return executable


def _run_synthetic_cli_workflow(cli: str, temp_dir: Path, *, env: dict[str, str]) -> None:
    fasta = temp_dir / "digest_input.fasta"
    fasta.write_text(">synthetic_digest\nAAAAGAATTCTTTT\n", encoding="utf-8")
    output = _run_output((cli, "digest", str(fasta), "--enzymes", "EcoRI"), env=env, cwd=temp_dir)
    if "EcoRI" not in output or "length" not in output:
        raise RuntimeError("synthetic digest workflow did not produce the expected CLI table")


def _venv_python(venv: Path) -> Path:
    if os.name == "nt":
        return venv / "Scripts" / "python.exe"
    return venv / "bin" / "python"


def _console_script(venv: Path, name: str) -> str:
    scripts = venv / ("Scripts" if os.name == "nt" else "bin")
    executable = scripts / (f"{name}.exe" if os.name == "nt" else name)
    if executable.exists():
        return str(executable)
    raise FileNotFoundError(f"Console script not found in smoke venv: {executable}")


def _isolated_python_env(extra: dict[str, str] | None = None) -> dict[str, str]:
    env = os.environ.copy()
    env.pop("PYTHONHOME", None)
    env.pop("PYTHONPATH", None)
    if extra:
        env.update(extra)
    return env


def _run(
    command: tuple[str, ...],
    *,
    env: dict[str, str] | None = None,
    cwd: Path | None = None,
    timeout: int | None = None,
) -> None:
    print("+ " + " ".join(command))
    subprocess.run(command, check=True, cwd=cwd, env=env, timeout=timeout)


def _run_output(
    command: tuple[str, ...],
    *,
    env: dict[str, str] | None = None,
    cwd: Path | None = None,
) -> str:
    print("+ " + " ".join(command))
    completed = subprocess.run(command, check=True, cwd=cwd, env=env, text=True, capture_output=True)
    output = completed.stdout + completed.stderr
    print(output)
    return output


if __name__ == "__main__":
    raise SystemExit(main())
