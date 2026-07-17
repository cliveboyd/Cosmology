param(
    [string]$RepoRoot = "C:\Users\clive\Documents\Cosmology",
    [string]$TaskRoot = "C:\Users\clive\Documents\Codex\2026-07-14\ok-continu"
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$outputs = Join-Path $TaskRoot "outputs"
$frPidFile = Join-Path $outputs "bbn_lithium_linx_chunked_resume.pid"
$frReport = Join-Path $RepoRoot "plamb_runs\diagnostics\bbn_lithium_linx_key_fr_network_overnight_20260716\bbn_lithium_linx_fr_network_interim_report.md"
$linxPython = Join-Path $RepoRoot "plamb_runs\envs\linx_py311\python.exe"
$bridgeProgram = Join-Path $RepoRoot "github_export\code\bbn_lithium\run_su2cmb_linx_thermal_bridge_2026-07-17.py"
$gateProgram = Join-Path $RepoRoot "github_export\code\bbn_lithium\analyze_su2_bbn_lithium_gate_2026-07-17.py"
$frontierProgram = Join-Path $RepoRoot "github_export\code\bbn_lithium\analyze_su2_bbn_lithium_frontier_2026-07-17.py"
$analysisPython = (Get-Command python.exe).Source

function Write-Stamped {
    param([string]$Message)
    Write-Output ("[{0}] {1}" -f (Get-Date -Format "o"), $Message)
}

foreach ($path in @($linxPython, $bridgeProgram, $gateProgram, $frontierProgram)) {
    if (-not (Test-Path -LiteralPath $path)) {
        throw "Required file not found: $path"
    }
}

Write-Stamped "SU2-CMB bridge queue started"
Write-Stamped "LINX Python: $linxPython"
Write-Stamped "Analysis Python: $analysisPython"

if (Test-Path -LiteralPath $frPidFile) {
    $frPidText = (Get-Content -LiteralPath $frPidFile -Raw).Trim()
    $frPid = 0
    if ([int]::TryParse($frPidText, [ref]$frPid)) {
        while ($true) {
            $process = Get-CimInstance Win32_Process -Filter "ProcessId=$frPid" -ErrorAction SilentlyContinue
            if ($null -eq $process) {
                break
            }
            if ($process.CommandLine -notlike "*resume_bbn_lithium_linx_key_chunks_2026-07-17.ps1*") {
                throw "FR PID $frPid was reused by an unrelated process; refusing to wait indefinitely."
            }
            Write-Stamped "Waiting for FR LINX wrapper PID $frPid to release memory"
            Start-Sleep -Seconds 30
        }
    }
}

Start-Sleep -Seconds 15
Write-Stamped "FR LINX wrapper is no longer active"

$frComplete = $false
if (Test-Path -LiteralPath $frReport) {
    $frText = Get-Content -LiteralPath $frReport -Raw
    $frComplete = (
        $frText -match "Successful unique points:\s*30503\s*/\s*30503" -and
        $frText -match "unresolved failed points:\s*0"
    )
}

if ($frComplete) {
    Write-Stamped "FR catalogue complete; regenerating registered gate"
    & $analysisPython -u $gateProgram
    if ($LASTEXITCODE -ne 0) {
        throw "Registered SU2 lithium gate failed with exit code $LASTEXITCODE"
    }
    Write-Stamped "Regenerating exploratory lithium frontier"
    & $analysisPython -u $frontierProgram
    if ($LASTEXITCODE -ne 0) {
        throw "SU2 lithium frontier failed with exit code $LASTEXITCODE"
    }
} else {
    Write-Stamped "WARNING: FR catalogue did not reach 30503/30503; bridge will proceed independently"
}

$env:XLA_PYTHON_CLIENT_PREALLOCATE = "false"
$env:MPLBACKEND = "Agg"
Write-Stamped "Starting literature-defined SU2-CMB LINX thermal bridge"
& $linxPython -u $bridgeProgram --network key --sampling-ntop 150 --rtol 1e-5 --atol 1e-9
if ($LASTEXITCODE -ne 0) {
    throw "SU2-CMB LINX bridge failed with exit code $LASTEXITCODE"
}

Write-Stamped "SU2-CMB LINX thermal bridge complete"
