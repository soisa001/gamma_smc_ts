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

$UvCommand = Get-Command uv -ErrorAction SilentlyContinue
if ($UvCommand) {
    $Uv = $UvCommand.Source
} else {
    Write-Host "Installing uv $UvVersion..."
    $Installer = Join-Path $Tools "uv-installer.ps1"
    Invoke-WebRequest `
        "https://releases.astral.sh/github/uv/releases/download/$UvVersion/uv-installer.ps1" `
        -OutFile $Installer
    & powershell -ExecutionPolicy Bypass -File $Installer
    $Uv = Join-Path $HOME ".local\bin\uv.exe"
    if (-not (Test-Path $Uv)) {
        $UvCommand = Get-Command uv -ErrorAction SilentlyContinue
        if (-not $UvCommand) { throw "uv installed but could not be located" }
        $Uv = $UvCommand.Source
    }
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
