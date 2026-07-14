$ErrorActionPreference = "Stop"

$TaskRoot = "C:\Users\clive\Documents\Codex\2026-07-14\ok-continu"
$CosmologyRoot = "C:\Users\clive\Documents\Cosmology"
$Python = "C:\Users\clive\anaconda3x\python.exe"
$Script = Join-Path $TaskRoot "work\su2_quaia_targeted_stability_mocks.py"
$OutDir = Join-Path $CosmologyRoot "plamb_runs\diagnostics\su2_quaia_targeted_stability_mocks_20260714"

if (-not (Test-Path -LiteralPath $Python)) { throw "Missing Python: $Python" }
if (-not (Test-Path -LiteralPath $Script)) { throw "Missing script: $Script" }

Push-Location $CosmologyRoot
try {
  & $Python $Script `
    --out $OutDir `
    --n-mocks 1500 `
    --progress-every 50 `
    --cache-randoms `
    --seed 290714
} finally {
  Pop-Location
}
