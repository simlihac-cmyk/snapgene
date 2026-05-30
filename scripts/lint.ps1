param(
    [string[]]$Paths = @("src", "tests"),
    [string]$Python = $env:PLASMIDLAB_PYTHON,
    [string]$VenvPath = ".venv-build",
    [switch]$ReuseVenv
)

$ErrorActionPreference = "Stop"
. (Join-Path $PSScriptRoot "common.ps1")

$hostPython = Resolve-PlasmidLabPython -Python $Python
$venvPython = Initialize-CleanVenv -Python $hostPython -VenvPath $VenvPath -Clean:(!$ReuseVenv)
Install-PlasmidLabDev -PythonPath $venvPython

$ruffArguments = @("-m", "ruff", "check") + $Paths
Invoke-CheckedCommand $venvPython $ruffArguments
