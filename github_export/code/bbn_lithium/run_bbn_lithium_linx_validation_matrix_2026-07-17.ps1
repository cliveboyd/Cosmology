param(
    [string]$RepoRoot = "C:\Users\clive\Documents\Cosmology",
    [string]$TaskOutputs = "C:\Users\clive\Documents\Codex\2026-07-14\ok-continu\outputs"
)

$ErrorActionPreference = "Stop"
$Python = Join-Path $RepoRoot "plamb_runs\envs\linx_py311\python.exe"
$Program = Join-Path $RepoRoot "github_export\code\bbn_lithium\validate_bbn_lithium_linx_matrix_2026-07-17.py"
$OutDir = Join-Path $RepoRoot "plamb_runs\diagnostics\bbn_lithium_linx_validation_matrix_20260717"
$Log = Join-Path $TaskOutputs "bbn_lithium_linx_validation_matrix_stdout.log"
$Err = Join-Path $TaskOutputs "bbn_lithium_linx_validation_matrix_stderr.log"

New-Item -ItemType Directory -Force -Path $OutDir, $TaskOutputs | Out-Null

function Write-RunLog([string]$Message) {
    $line = "[$(Get-Date -Format o)] $Message"
    $line | Tee-Object -FilePath $Log -Append
}

try {
    $matrixPath = Join-Path $OutDir "validation_matrix.csv"
    $configPath = Join-Path $OutDir "validation_matrix_config.json"
    $registrationPath = Join-Path $OutDir "validation_matrix_registration.md"
    if ((Test-Path -LiteralPath $matrixPath) -and (Test-Path -LiteralPath $configPath) -and (Test-Path -LiteralPath $registrationPath)) {
        Write-RunLog "Using the existing frozen preregistered validation matrix"
    }
    else {
        Write-RunLog "Preparing preregistered validation matrix"
        & $Python $Program --mode prepare --outdir $OutDir 2>> $Err | Tee-Object -FilePath $Log -Append
        if ($LASTEXITCODE -ne 0) { throw "Preparation failed with exit code $LASTEXITCODE" }
    }

    foreach ($network in @("key", "small", "full")) {
        Write-RunLog "Starting $network network"
        & $Python $Program --mode run --network $network --outdir $OutDir --resume 2>> $Err | Tee-Object -FilePath $Log -Append
        if ($LASTEXITCODE -ne 0) { throw "$network network failed with exit code $LASTEXITCODE" }
        Write-RunLog "Completed $network network"
    }

    Write-RunLog "Analysing completed matrix"
    & $Python $Program --mode analyse --outdir $OutDir 2>> $Err | Tee-Object -FilePath $Log -Append
    if ($LASTEXITCODE -ne 0) { throw "Analysis failed with exit code $LASTEXITCODE" }
    Write-RunLog "Validation matrix complete"
}
catch {
    Write-RunLog "FAILED: $($_.Exception.Message)"
    throw
}
