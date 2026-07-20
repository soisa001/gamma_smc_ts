$ErrorActionPreference = "Stop"
$Repo = Split-Path -Parent $PSScriptRoot
$UvCommand = Get-Command uv -ErrorAction SilentlyContinue
$Uv = if ($UvCommand) { $UvCommand.Source } else { Join-Path $Repo ".tools\uv-bin\uv.exe" }
if (-not (Test-Path $Uv) -and -not $UvCommand) {
    throw "uv is not installed; run scripts\bootstrap_uv.ps1 first."
}

$LocalSlim = Join-Path $Repo ".native\Library\bin\slim.exe"
if (-not $env:SLIM_BIN -and (Test-Path $LocalSlim)) { $env:SLIM_BIN = $LocalSlim }
$LocalDecoder = Join-Path $Repo "bin\gamma_smc"
if (-not $env:GAMMA_SMC_BIN -and (Test-Path $LocalDecoder)) {
    $env:GAMMA_SMC_BIN = $LocalDecoder
}

& $Uv run --project $Repo --frozen --all-extras gamma-smc-aou @args
exit $LASTEXITCODE
