"""PyInstaller entrypoint with headless release-smoke commands."""

from __future__ import annotations

import json
import sys
from importlib.util import find_spec
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as package_version
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    """Handle smoke flags before importing the full Qt main window."""

    arguments = list(sys.argv[1:] if argv is None else argv)
    if "--help" in arguments or "-h" in arguments:
        _write_stdout(_gui_cli_help())
        return 0
    if "--version" in arguments:
        _write_stdout(f"PlasmidLab {_plasmidlab_version()}\n")
        return 0
    if "--smoke-json" in arguments:
        index = arguments.index("--smoke-json")
        try:
            output_path = Path(arguments[index + 1])
        except IndexError as error:
            msg = "--smoke-json requires an output path"
            raise SystemExit(msg) from error

        resource_entries, feature_library_entries = _bundled_feature_library_counts()
        package_resources_available = resource_entries > 0
        feature_library_available = feature_library_entries > 0
        output_path.write_text(
            json.dumps(
                {
                    "application": "PlasmidLab",
                    "cli": True,
                    "feature_library": feature_library_available,
                    "version": _plasmidlab_version(),
                    "package_resources": package_resources_available,
                    "feature_library_entries": feature_library_entries,
                    "feature_resource_entries": resource_entries,
                    "ok": package_resources_available and feature_library_available,
                    "python_version_ok": sys.version_info[:2] >= (3, 12),
                    "qt_import": find_spec("PySide6") is not None,
                },
                sort_keys=True,
            ),
            encoding="utf-8",
        )
        return 0

    from plasmidlab.gui.main_window import main as gui_main

    return gui_main(arguments)


def _plasmidlab_version() -> str:
    try:
        return package_version("plasmidlab")
    except PackageNotFoundError:
        return "0.1.0"


def _bundled_feature_library_counts() -> tuple[int, int]:
    resource_count = 0
    feature_count = 0
    for directory in _candidate_feature_library_dirs():
        if not directory.is_dir():
            continue
        for resource in sorted(directory.glob("*.json"), key=lambda path: path.name):
            resource_count += 1
            data = json.loads(resource.read_text(encoding="utf-8"))
            if isinstance(data, dict):
                raw_entries = data.get("features", [data])
            elif isinstance(data, list):
                raw_entries = data
            else:
                raw_entries = []
            feature_count += sum(
                1
                for entry in raw_entries
                if isinstance(entry, dict) and entry.get("name") and entry.get("sequence")
            )
        if resource_count:
            break
    return resource_count, feature_count


def _candidate_feature_library_dirs() -> tuple[Path, ...]:
    candidates: list[Path] = []
    frozen_root = getattr(sys, "_MEIPASS", None)
    if frozen_root:
        candidates.append(Path(frozen_root) / "plasmidlab" / "data" / "features")
    if getattr(sys, "frozen", False):
        executable_dir = Path(sys.executable).resolve().parent
        candidates.extend(
            (
                executable_dir / "_internal" / "plasmidlab" / "data" / "features",
                executable_dir / "plasmidlab" / "data" / "features",
            )
        )
    candidates.append(Path(__file__).resolve().parents[1] / "data" / "features")
    return tuple(dict.fromkeys(candidates))


def _gui_cli_help() -> str:
    return (
        "usage: PlasmidLab [--help] [--version] [--smoke-json PATH]\n\n"
        "Open the PlasmidLab desktop GUI.\n\n"
        "options:\n"
        "  --help             show this help and exit\n"
        "  --version          print the PlasmidLab version and exit\n"
        "  --smoke-json PATH  write a headless frozen-app smoke report and exit\n"
    )


def _write_stdout(text: str) -> None:
    stream = getattr(sys, "stdout", None)
    if stream is not None:
        stream.write(text)


if __name__ == "__main__":
    raise SystemExit(main())
