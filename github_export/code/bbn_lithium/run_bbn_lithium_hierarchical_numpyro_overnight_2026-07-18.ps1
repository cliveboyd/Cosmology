$ErrorActionPreference = "Stop"

throw "Disabled after preregistered smoke tests: reverse-mode LINX gradients produced a non-finite implicit solve and forward-mode is unsupported by LINX custom_vjp. Use run_bbn_lithium_hierarchical_surrogate_pipeline_2026-07-18.ps1."

$python = "C:\Users\clive\Documents\Cosmology\plamb_runs\envs\linx_py311\python.exe"
$script = "C:\Users\clive\Documents\Cosmology\github_export\code\bbn_lithium\run_bbn_lithium_hierarchical_numpyro_2026-07-18.py"

Write-Output "[$(Get-Date -Format s)] Starting preregistered hierarchical LINX BBN posterior"
& $python $script --mode overnight --arm both
if ($LASTEXITCODE -ne 0) {
    throw "Hierarchical LINX BBN posterior failed with exit code $LASTEXITCODE"
}
Write-Output "[$(Get-Date -Format s)] Hierarchical LINX BBN posterior and analysis complete"
