"""PyInstaller entrypoint with headless release-smoke commands."""

from __future__ import annotations

import json
import sys
from importlib import resources
from importlib.util import find_spec
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    """Handle smoke flags before importing the full Qt main window."""

    arguments = list(sys.argv[1:] if argv is None else argv)
    if "--help" in arguments or "-h" in arguments:
        _write_stdout(_gui_cli_help())
        return 0
    if "--version" in arguments:
        import plasmidlab

        _write_stdout(f"PlasmidLab {plasmidlab.__version__}\n")
        return 0
    if "--smoke-json" in arguments:
        index = arguments.index("--smoke-json")
        try:
            output_path = Path(arguments[index + 1])
        except IndexError as error:
            msg = "--smoke-json requires an output path"
            raise SystemExit(msg) from error

        import plasmidlab
        from plasmidlab.core.feature_annotation import load_feature_library

        feature_resource_root = resources.files("plasmidlab.data.features")
        resource_entries = [
            resource for resource in feature_resource_root.iterdir() if resource.name.endswith(".json")
        ]
        library = load_feature_library()
        package_resources_available = bool(resource_entries)
        feature_library_available = bool(library)
        output_path.write_text(
            json.dumps(
                {
                    "application": "PlasmidLab",
                    "cli": True,
                    "feature_library": feature_library_available,
                    "version": plasmidlab.__version__,
                    "package_resources": package_resources_available,
                    "feature_library_entries": len(library),
                    "feature_resource_entries": len(resource_entries),
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
