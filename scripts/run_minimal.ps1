param(
    [int]$Replicates = 200,
    [int]$Diploids = 100,
    [int]$Length = 200000,
    [int]$Workers = 20,
    [string]$OutputDir = "sim_results/minimal"
)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
$env:PYTHONPATH = Join-Path $repo "python"
$python = if (Test-Path (Join-Path (Split-Path $repo -Parent) ".venv/Scripts/python.exe")) {
    Join-Path (Split-Path $repo -Parent) ".venv/Scripts/python.exe"
} else { "python" }

& $python -m gamma_smc_aou.cli simulate `
    --output-dir (Join-Path $repo $OutputDir) `
    --replicates $Replicates --diploids $Diploids --length $Length `
    --ne 10000 --mutation-rate 1.25e-8 --recombination-rate 1e-8 --seed 1729 --workers $Workers

& $python -m gamma_smc_aou.cli validate-null `
    --sim-glob "$repo/$OutputDir/truth_summaries/*.tsv" `
    --sequence-length $Length --output-dir "$repo/$OutputDir/calibration"

Write-Host "Minimal fixed-standard-coalescent results: $repo/$OutputDir"
