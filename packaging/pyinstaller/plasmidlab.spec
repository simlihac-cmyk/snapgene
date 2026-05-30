# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller build configuration for the PlasmidLab desktop app."""

from __future__ import annotations

import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_data_files, collect_submodules


ROOT = Path(SPECPATH).parents[1]
SRC = ROOT / "src"
ICON = ROOT / "assets" / "icons" / "plasmidlab.ico"

sys.path.insert(0, str(SRC))

datas = [
    (str(ROOT / "assets" / "icons" / "plasmidlab.svg"), "assets/icons"),
]
datas += collect_data_files("plasmidlab.data")
datas += collect_data_files("plasmidlab.resources")

hiddenimports = []
hiddenimports += [
    "Bio.Align",
    "Bio.Data.CodonTable",
    "Bio.Seq",
    "Bio.SeqFeature",
    "Bio.SeqIO",
    "Bio.SeqIO.AbiIO",
    "Bio.SeqIO.FastaIO",
    "Bio.SeqIO.InsdcIO",
    "Bio.SeqRecord",
]
hiddenimports += collect_submodules("Bio.Restriction")
hiddenimports += collect_submodules("primer3")

a = Analysis(
    [str(SRC / "plasmidlab" / "gui" / "main_window.py")],
    pathex=[str(SRC)],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="PlasmidLab",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=str(ICON) if ICON.exists() else None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name="PlasmidLab",
)
