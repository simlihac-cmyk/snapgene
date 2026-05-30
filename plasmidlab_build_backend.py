"""PEP 517 backend guard for PlasmidLab release builds.

The project metadata requires Python 3.12+, but some build frontends can still create
metadata under an older interpreter when a maintainer bypasses the release scripts. This
backend wrapper fails before delegating to Hatchling so bare ``python -m build`` obeys
the same version floor.
"""

from __future__ import annotations

import sys
from collections.abc import Mapping
from typing import Any


MIN_VERSION = (3, 12)
MIN_VERSION_LABEL = "3.12"


def _ensure_supported_python(version_info: tuple[int, ...] | None = None) -> None:
    version = tuple(version_info or sys.version_info[:3])
    if version[:2] >= MIN_VERSION:
        return
    detected = ".".join(str(part) for part in version[:3])
    msg = (
        f"PlasmidLab requires Python {MIN_VERSION_LABEL} or newer for builds; "
        f"detected Python {detected}. Create a clean Python {MIN_VERSION_LABEL}+ "
        "virtual environment or set PLASMIDLAB_PYTHON to a supported interpreter."
    )
    raise RuntimeError(msg)


def _hatchling_build() -> Any:
    from hatchling import build

    return build


def get_requires_for_build_wheel(config_settings: Mapping[str, Any] | None = None) -> list[str]:
    _ensure_supported_python()
    return _hatchling_build().get_requires_for_build_wheel(config_settings)


def get_requires_for_build_sdist(config_settings: Mapping[str, Any] | None = None) -> list[str]:
    _ensure_supported_python()
    return _hatchling_build().get_requires_for_build_sdist(config_settings)


def prepare_metadata_for_build_wheel(
    metadata_directory: str,
    config_settings: Mapping[str, Any] | None = None,
) -> str:
    _ensure_supported_python()
    return _hatchling_build().prepare_metadata_for_build_wheel(metadata_directory, config_settings)


def build_wheel(
    wheel_directory: str,
    config_settings: Mapping[str, Any] | None = None,
    metadata_directory: str | None = None,
) -> str:
    _ensure_supported_python()
    return _hatchling_build().build_wheel(wheel_directory, config_settings, metadata_directory)


def build_sdist(
    sdist_directory: str,
    config_settings: Mapping[str, Any] | None = None,
) -> str:
    _ensure_supported_python()
    return _hatchling_build().build_sdist(sdist_directory, config_settings)


def get_requires_for_build_editable(config_settings: Mapping[str, Any] | None = None) -> list[str]:
    _ensure_supported_python()
    return _hatchling_build().get_requires_for_build_editable(config_settings)


def prepare_metadata_for_build_editable(
    metadata_directory: str,
    config_settings: Mapping[str, Any] | None = None,
) -> str:
    _ensure_supported_python()
    return _hatchling_build().prepare_metadata_for_build_editable(metadata_directory, config_settings)


def build_editable(
    wheel_directory: str,
    config_settings: Mapping[str, Any] | None = None,
    metadata_directory: str | None = None,
) -> str:
    _ensure_supported_python()
    return _hatchling_build().build_editable(wheel_directory, config_settings, metadata_directory)
