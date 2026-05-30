param(
    [string]$Python = $env:PLASMIDLAB_PYTHON,
    [string]$VenvPath = ".venv-build",
    [switch]$ReuseVenv
)

$ErrorActionPreference = "Stop"
. (Join-Path $PSScriptRoot "common.ps1")

$hostPython = Resolve-PlasmidLabPython -Python $Python
$venvPython = Initialize-CleanVenv -Python $hostPython -VenvPath $VenvPath -Clean:(!$ReuseVenv)
Install-PlasmidLabDev -PythonPath $venvPython

Invoke-CheckedCommand $venvPython @("scripts/generate_icon.py")
Invoke-CheckedCommand $venvPython @(
    "-m",
    "PyInstaller",
    "--noconfirm",
    "build_support/pyinstaller/plasmidlab.spec"
)
