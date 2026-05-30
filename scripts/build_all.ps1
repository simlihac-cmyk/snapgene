param(
    [string]$Python = $env:PLASMIDLAB_PYTHON,
    [string]$VenvPath = ".venv-build"
)

$ErrorActionPreference = "Stop"
. (Join-Path $PSScriptRoot "common.ps1")

$hostPython = Resolve-PlasmidLabPython -Python $Python
Remove-ReleaseArtifacts
$venvPython = Initialize-CleanVenv -Python $hostPython -VenvPath $VenvPath -Clean
Install-PlasmidLabDev -PythonPath $venvPython

Invoke-CheckedCommand $venvPython @("-m", "pytest")
Invoke-CheckedCommand $venvPython @("-m", "ruff", "check", "src", "tests", "scripts")
Invoke-CheckedCommand $venvPython @("scripts/generate_icon.py")
Invoke-CheckedCommand $venvPython @("-m", "build", "--no-isolation")
Invoke-CheckedCommand $venvPython @("scripts/smoke_release.py", "--wheel-only")
Invoke-CheckedCommand $venvPython @(
    "-m",
    "PyInstaller",
    "--noconfirm",
    "packaging/pyinstaller/plasmidlab.spec"
)
Invoke-CheckedCommand $venvPython @("scripts/smoke_release.py", "--pyinstaller-only")
Invoke-CheckedCommand $venvPython @("scripts/generate_checksums.py")
