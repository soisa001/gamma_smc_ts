[CmdletBinding()]
param(
    [switch]$Full,
    [switch]$SkipNative,
    [switch]$SkipTests,
    [string]$PythonVersion = "3.11"
)

$ErrorActionPreference = "Stop"
$Repo = Split-Path -Parent $PSScriptRoot
$Tools = Join-Path $Repo ".tools"
$Native = Join-Path $Repo ".native"
$UvVersion = "0.11.16"
$MicromambaVersion = "2.6.2-1"
New-Item -ItemType Directory -Force -Path $Tools | Out-Null
$RunningWindows = $env:OS -eq "Windows_NT"

function Invoke-Checked {
    param([string]$File, [string[]]$Arguments)
    & $File @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Command failed ($LASTEXITCODE): $File $($Arguments -join ' ')"
    }
}

$LocalUv = Join-Path $Tools "uv-bin\uv.exe"
$UvCommand = Get-Command uv -ErrorAction SilentlyContinue
$Uv = if ($UvCommand) { $UvCommand.Source } else { $null }
$DetectedUvVersion = if ($Uv) { (& $Uv --version 2>$null) } else { $null }
if ($DetectedUvVersion -ne "uv $UvVersion") {
    if ($Uv) {
        Write-Host "Ignoring incompatible $DetectedUvVersion at $Uv."
    }
    if ((Test-Path $LocalUv) -and ((& $LocalUv --version 2>$null) -eq "uv $UvVersion")) {
        $Uv = $LocalUv
    } else {
        Write-Host "Installing repository-pinned uv $UvVersion..."
        New-Item -ItemType Directory -Force -Path (Split-Path -Parent $LocalUv) | Out-Null
        $env:UV_INSTALL_DIR = Split-Path -Parent $LocalUv
        $Installer = Join-Path $Tools "uv-installer.ps1"
        Invoke-WebRequest `
            "https://releases.astral.sh/github/uv/releases/download/$UvVersion/uv-installer.ps1" `
            -OutFile $Installer
        & powershell -ExecutionPolicy Bypass -File $Installer -NoModifyPath
        $Uv = $LocalUv
    }
}
if ((& $Uv --version 2>$null) -ne "uv $UvVersion") {
    throw "Failed to select uv $UvVersion at $Uv"
}

if ($Full) {
    throw "The Gamma-SMC decoder is Linux x86_64/AVX2 software. Run this repository under Linux/WSL for -Full; omit -Full for simulations on Windows/macOS."
}

$Slim = $env:SLIM_BIN
if (-not $SkipNative) {
    if (-not $RunningWindows) {
        throw "Use scripts/bootstrap_uv.sh on Linux or macOS."
    }
    $Mamba = Join-Path $Tools "micromamba.exe"
    if (-not (Test-Path $Mamba)) {
        Write-Host "Installing micromamba $MicromambaVersion..."
        Invoke-WebRequest `
            "https://github.com/mamba-org/micromamba-releases/releases/download/$MicromambaVersion/micromamba-win-64" `
            -OutFile $Mamba
    }
    $MambaArgs = @("-y", "-p", $Native, "-c", "conda-forge", "--strict-channel-priority", "slim=5.2")
    if (Test-Path (Join-Path $Native "conda-meta")) {
        Invoke-Checked $Mamba (@("install") + $MambaArgs)
    } else {
        Invoke-Checked $Mamba (@("create") + $MambaArgs)
    }
    $Slim = Join-Path $Native "Library\bin\slim.exe"
}

Write-Host "Installing locked Python environment..."
Invoke-Checked $Uv @("python", "install", $PythonVersion)
Invoke-Checked $Uv @("sync", "--project", $Repo, "--frozen", "--all-extras", "--python", $PythonVersion)

if ($Slim) { $env:SLIM_BIN = $Slim }
$VerifyArgs = @("run", "--project", $Repo, "--frozen", "--all-extras", "python", (Join-Path $PSScriptRoot "verify_install.py"), "--require-slim")
if ($Full) { $VerifyArgs += "--require-decoder" }
Invoke-Checked $Uv $VerifyArgs

if (-not $SkipTests) {
    Write-Host "Running test suite..."
    Invoke-Checked $Uv @("run", "--project", $Repo, "--frozen", "--all-extras", "pytest", "-q")
}

Write-Host ""
Write-Host "Ready. Run commands through: scripts\aou.ps1 <gamma-smc-aou arguments>"
