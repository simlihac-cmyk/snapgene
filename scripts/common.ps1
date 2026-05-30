$ErrorActionPreference = "Stop"

$script:RepoRoot = Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $script:RepoRoot

function Invoke-CheckedCommand {
    param(
        [Parameter(Mandatory = $true)][string]$Command,
        [string[]]$Arguments = @()
    )

    & $Command @Arguments
    if ($LASTEXITCODE -ne 0) {
        exit $LASTEXITCODE
    }
}

function Invoke-PythonCommand {
    param(
        [Parameter(Mandatory = $true)]$Python,
        [string[]]$Arguments = @()
    )

    & $Python.FilePath @($Python.PrefixArgs + $Arguments)
    if ($LASTEXITCODE -ne 0) {
        exit $LASTEXITCODE
    }
}

function Test-PythonCandidate {
    param(
        [Parameter(Mandatory = $true)][string]$FilePath,
        [string[]]$PrefixArgs = @()
    )

    $checkScript = Join-Path $script:RepoRoot "scripts/check_python.py"
    try {
        & $FilePath @($PrefixArgs + @($checkScript)) 1>$null 2>$null
        return $LASTEXITCODE -eq 0
    }
    catch {
        return $false
    }
}

function Resolve-PlasmidLabPython {
    param(
        [string]$Python = $env:PLASMIDLAB_PYTHON
    )

    $candidates = @()
    if ($Python) {
        $candidates += [pscustomobject]@{ FilePath = $Python; PrefixArgs = @() }
    }
    $candidates += [pscustomobject]@{ FilePath = "python"; PrefixArgs = @() }
    $candidates += [pscustomobject]@{ FilePath = "python3.13"; PrefixArgs = @() }
    $candidates += [pscustomobject]@{ FilePath = "python3.12"; PrefixArgs = @() }
    $candidates += [pscustomobject]@{ FilePath = "py"; PrefixArgs = @("-3.13") }
    $candidates += [pscustomobject]@{ FilePath = "py"; PrefixArgs = @("-3.12") }

    foreach ($candidate in $candidates) {
        if (-not (Get-Command $candidate.FilePath -ErrorAction SilentlyContinue)) {
            continue
        }
        if (Test-PythonCandidate -FilePath $candidate.FilePath -PrefixArgs $candidate.PrefixArgs) {
            return $candidate
        }
    }

    if (Get-Command "python" -ErrorAction SilentlyContinue) {
        & python (Join-Path $script:RepoRoot "scripts/check_python.py")
    }

    $message = @(
        "No supported Python interpreter found.",
        "PlasmidLab release tooling requires Python 3.12 or newer.",
        "Install Python 3.12+ and rerun, or set PLASMIDLAB_PYTHON to its python.exe path."
    ) -join " "
    Write-Host $message -ForegroundColor Red
    exit 1
}

function Get-VenvPythonPath {
    param(
        [Parameter(Mandatory = $true)][string]$VenvPath
    )

    if ($IsWindows -or $env:OS -eq "Windows_NT") {
        return Join-Path $VenvPath "Scripts/python.exe"
    }
    return Join-Path $VenvPath "bin/python"
}

function Assert-SafeVenvPath {
    param(
        [Parameter(Mandatory = $true)][string]$VenvPath
    )

    $root = [System.IO.Path]::GetFullPath($script:RepoRoot)
    $target = [System.IO.Path]::GetFullPath((Join-Path $script:RepoRoot $VenvPath))
    $rootWithSeparator = $root.TrimEnd([System.IO.Path]::DirectorySeparatorChar) + [System.IO.Path]::DirectorySeparatorChar
    $leaf = Split-Path -Leaf $target
    if (-not $target.StartsWith($rootWithSeparator, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to manage virtual environment outside repository: $target"
    }
    if (-not $leaf.StartsWith(".venv", [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to clean a virtual environment path that does not start with .venv: $target"
    }
    return $target
}

function Remove-ReleaseArtifacts {
    $allowedLeaves = @("build", "dist")
    $root = [System.IO.Path]::GetFullPath($script:RepoRoot)
    $rootWithSeparator = $root.TrimEnd([System.IO.Path]::DirectorySeparatorChar) + [System.IO.Path]::DirectorySeparatorChar
    foreach ($leaf in $allowedLeaves) {
        $target = [System.IO.Path]::GetFullPath((Join-Path $script:RepoRoot $leaf))
        if (-not $target.StartsWith($rootWithSeparator, [System.StringComparison]::OrdinalIgnoreCase)) {
            throw "Refusing to remove artifact path outside repository: $target"
        }
        if ((Split-Path -Leaf $target) -notin $allowedLeaves) {
            throw "Refusing to remove unexpected artifact path: $target"
        }
        if (Test-Path -LiteralPath $target) {
            Remove-Item -LiteralPath $target -Recurse -Force
        }
    }
}

function Initialize-CleanVenv {
    param(
        [Parameter(Mandatory = $true)]$Python,
        [string]$VenvPath = ".venv-build",
        [switch]$Clean
    )

    $absoluteVenv = Assert-SafeVenvPath -VenvPath $VenvPath
    if ($Clean -and (Test-Path -LiteralPath $absoluteVenv)) {
        Remove-Item -LiteralPath $absoluteVenv -Recurse -Force
    }
    if (-not (Test-Path -LiteralPath $absoluteVenv)) {
        Invoke-PythonCommand $Python @("-m", "venv", $absoluteVenv)
    }

    $venvPython = Get-VenvPythonPath -VenvPath $absoluteVenv
    Invoke-CheckedCommand $venvPython @((Join-Path $script:RepoRoot "scripts/check_python.py"))
    Invoke-CheckedCommand $venvPython @("-m", "pip", "install", "--upgrade", "pip")
    return $venvPython
}

function Install-PlasmidLabDev {
    param(
        [Parameter(Mandatory = $true)][string]$PythonPath
    )

    $constraints = Join-Path $script:RepoRoot "constraints.txt"
    $editableSpec = "$script:RepoRoot[dev]"
    Invoke-CheckedCommand $PythonPath @(
        "-m",
        "pip",
        "install",
        "--constraint",
        $constraints,
        "-e",
        $editableSpec
    )
}
