"""Check that PlasmidLab tooling is running on a supported Python version."""

from __future__ import annotations

import sys
from collections.abc import Sequence


MIN_VERSION = (3, 12)
MIN_VERSION_LABEL = "3.12"


def version_label(version_info: Sequence[int]) -> str:
    """Return a short major.minor.patch label for a version tuple."""

    major, minor, micro = version_info[:3]
    return f"{major}.{minor}.{micro}"


def is_supported(version_info: Sequence[int]) -> bool:
    return tuple(version_info[:2]) >= MIN_VERSION


def main() -> int:
    detected = version_label(sys.version_info)
    executable = sys.executable or "<unknown executable>"
    print(f"PlasmidLab Python check: detected Python {detected} at {executable}")
    if not is_supported(sys.version_info):
        print(
            f"ERROR: PlasmidLab requires Python {MIN_VERSION_LABEL} or newer. "
            f"Create a clean Python {MIN_VERSION_LABEL}+ virtual environment and rerun this command.",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
