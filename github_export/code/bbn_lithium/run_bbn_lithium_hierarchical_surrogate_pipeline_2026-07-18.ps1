$ErrorActionPreference = "Stop"

$linxPython = "C:\Users\clive\Documents\Cosmology\plamb_runs\envs\linx_py311\python.exe"
$analysisPython = "C:\Users\clive\anaconda3x\python.exe"
$exactProgramme = "C:\Users\clive\Documents\Cosmology\github_export\code\bbn_lithium\run_bbn_lithium_hierarchical_exact_design_2026-07-18.py"
$analysisProgramme = "C:\Users\clive\Documents\Cosmology\github_export\code\bbn_lithium\analyse_bbn_lithium_hierarchical_surrogate_2026-07-18.py"
$runDir = "C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\bbn_lithium_hierarchical_surrogate_20260718"
$validationPoints = Join-Path $runDir "surrogate_posterior_exact_validation_points.csv"
$validationResults = Join-Path $runDir "exact_linx_posterior_validation_results.csv"

function Invoke-Checked {
    param([string]$Label, [scriptblock]$Command)
    Write-Output "[$(Get-Date -Format o)] START $Label"
    & $Command
    if ($LASTEXITCODE -ne 0) {
        throw "$Label failed with exit code $LASTEXITCODE"
    }
    Write-Output "[$(Get-Date -Format o)] DONE $Label"
}

Invoke-Checked "exact LINX 4800-point design" {
    & $linxPython $exactProgramme --mode design --summary-every 25
}

Invoke-Checked "surrogate training and gradient-free posterior" {
    & $analysisPython $analysisProgramme --mode train-sample --walkers 64 --burn 3000 --steps 12000 --validation-per-arm 160
}

Invoke-Checked "exact LINX posterior-draw validation" {
    & $linxPython $exactProgramme --mode validate --input-points $validationPoints --output-csv $validationResults --summary-every 20
}

Invoke-Checked "final exact-corrected gate analysis" {
    & $analysisPython $analysisProgramme --mode finalise --validation-per-arm 160
}

Write-Output "[$(Get-Date -Format o)] Hierarchical surrogate pipeline complete"
