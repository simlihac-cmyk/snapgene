from __future__ import annotations

import sys

import pytest

from plasmidlab.gui.main_window import main as gui_main
from plasmidlab.gui.frozen_entry import main as frozen_main
from scripts import generate_checksums
from scripts.smoke_release import _isolated_python_env


def test_smoke_release_isolated_env_removes_source_import_paths(monkeypatch) -> None:
    monkeypatch.setenv("PYTHONHOME", "should-not-leak")
    monkeypatch.setenv("PYTHONPATH", "src")

    env = _isolated_python_env()

    assert "PYTHONHOME" not in env
    assert "PYTHONPATH" not in env


def test_gui_smoke_json_reports_package_resource_status(tmp_path) -> None:
    smoke_path = tmp_path / "smoke.json"

    assert gui_main(["--smoke-json", str(smoke_path)]) == 0

    smoke_json = smoke_path.read_text(encoding="utf-8")
    assert '"feature_library": true' in smoke_json
    assert '"package_resources": true' in smoke_json
    assert '"feature_library_entries":' in smoke_json


def test_frozen_entry_smoke_json_reports_package_resource_status(tmp_path) -> None:
    smoke_path = tmp_path / "frozen-smoke.json"

    assert frozen_main(["--smoke-json", str(smoke_path)]) == 0

    smoke_json = smoke_path.read_text(encoding="utf-8")
    assert '"feature_library": true' in smoke_json
    assert '"package_resources": true' in smoke_json
    assert '"feature_library_entries":' in smoke_json


def test_generate_checksums_requires_all_release_artifacts(tmp_path, monkeypatch) -> None:
    dist = tmp_path / "dist"
    pyinstaller_dir = dist / "PlasmidLab"
    pyinstaller_dir.mkdir(parents=True)

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "generate_checksums.py",
            "--dist",
            str(dist),
            "--pyinstaller-dir",
            str(pyinstaller_dir),
        ],
    )

    with pytest.raises(FileNotFoundError, match=r"wheel \(\*\.whl\).*sdist"):
        generate_checksums.main()


def test_generate_checksums_writes_checksums_for_required_artifacts(tmp_path, monkeypatch) -> None:
    dist = tmp_path / "dist"
    pyinstaller_dir = dist / "PlasmidLab"
    pyinstaller_dir.mkdir(parents=True)
    (dist / "plasmidlab-0.1.0-py3-none-any.whl").write_text("wheel", encoding="utf-8")
    (dist / "plasmidlab-0.1.0.tar.gz").write_text("sdist", encoding="utf-8")
    (pyinstaller_dir / "PlasmidLab.exe").write_text("frozen app", encoding="utf-8")

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "generate_checksums.py",
            "--dist",
            str(dist),
            "--pyinstaller-dir",
            str(pyinstaller_dir),
        ],
    )

    assert generate_checksums.main() == 0
    checksums = (dist / "SHA256SUMS.txt").read_text(encoding="utf-8")
    assert "plasmidlab-0.1.0-py3-none-any.whl" in checksums
    assert "plasmidlab-0.1.0.tar.gz" in checksums
    assert "PlasmidLab-windows.zip" in checksums
