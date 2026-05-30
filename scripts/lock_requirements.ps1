param(
    [string]$Python = $env:PLASMIDLAB_PYTHON,
    [string]$VenvPath = ".venv-lock"
)

$ErrorActionPreference = "Stop"
. (Join-Path $PSScriptRoot "common.ps1")

$hostPython = Resolve-PlasmidLabPython -Python $Python
$venvPython = Initialize-CleanVenv -Python $hostPython -VenvPath $VenvPath -Clean
$editableSpec = "$script:RepoRoot[dev]"
Invoke-CheckedCommand $venvPython @("-m", "pip", "install", "-e", $editableSpec)
$lockPath = Join-Path $script:RepoRoot "requirements.lock"
& $venvPython -m pip freeze --all | Set-Content -LiteralPath $lockPath -Encoding UTF8
if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
}

$message = @(
    "Wrote requirements.lock from a clean Python 3.12+ environment.",
    "Review it and copy release-build pins into constraints.txt before tagging."
) -join " "
Write-Host $message
