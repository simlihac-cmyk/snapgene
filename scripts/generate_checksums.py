"""Generate SHA256 checksums for PlasmidLab release artifacts."""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dist", type=Path, default=Path("dist"))
    parser.add_argument("--pyinstaller-dir", type=Path, default=Path("dist/PlasmidLab"))
    args = parser.parse_args()

    args.dist.mkdir(parents=True, exist_ok=True)
    wheels = sorted(args.dist.glob("*.whl"))
    sdists = sorted(args.dist.glob("*.tar.gz"))
    missing = []
    if not wheels:
        missing.append("wheel (*.whl)")
    if not sdists:
        missing.append("sdist (*.tar.gz)")
    if not args.pyinstaller_dir.exists():
        missing.append(f"PyInstaller directory ({args.pyinstaller_dir})")
    if missing:
        joined = ", ".join(missing)
        raise FileNotFoundError(f"Missing required release artifact(s): {joined}")

    pyinstaller_zip = _zip_pyinstaller_artifact(args.pyinstaller_dir)
    artifacts = [*wheels, *sdists, pyinstaller_zip]

    checksum_path = args.dist / "SHA256SUMS.txt"
    lines = [f"{_sha256(path)}  {path.name}" for path in artifacts]
    checksum_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {checksum_path}")
    return 0


def _zip_pyinstaller_artifact(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"PyInstaller artifact directory not found: {path}")
    platform_label = "windows" if os.name == "nt" else "desktop"
    archive_base = path.parent / f"PlasmidLab-{platform_label}"
    archive_path = archive_base.with_suffix(".zip")
    if archive_path.exists():
        archive_path.unlink()
    created = shutil.make_archive(str(archive_base), "zip", root_dir=path.parent, base_dir=path.name)
    return Path(created)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


if __name__ == "__main__":
    raise SystemExit(main())
