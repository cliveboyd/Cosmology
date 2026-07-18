$ErrorActionPreference = "Stop"

$root = "C:\Users\clive\Documents\Cosmology"
$pipeline = Join-Path $root "github_export\code\quaia_su2\run_sdss_dr16q_v4_cross_catalogue_pipeline_2026-07-18.py"
$dataDir = Join-Path $root "external_datasets\sdss_dr16q_v4"
$pidFile = Join-Path $dataDir "sdss_dr16q_v4_cross_catalogue_pipeline.pid"
$completedReport = Join-Path $root "github_export\results\2026-07-18\su2_quaia\sdss_dr16q_v4_cross_catalogue_validation_2026-07-18_report.md"
$logDir = "C:\Users\clive\Documents\Codex\2026-07-14\ok-continu\outputs"
$stdoutLog = Join-Path $logDir "sdss_dr16q_v4_cross_catalogue_pipeline_stdout.log"
$stderrLog = Join-Path $logDir "sdss_dr16q_v4_cross_catalogue_pipeline_stderr.log"

New-Item -ItemType Directory -Force -Path $dataDir | Out-Null
New-Item -ItemType Directory -Force -Path $logDir | Out-Null

if (Test-Path -LiteralPath $completedReport) {
    Write-Host "Validation report already exists; no duplicate pipeline was started."
    Write-Host "Report: $completedReport"
    exit 0
}

if (Test-Path -LiteralPath $pidFile) {
    $existingPidText = (Get-Content -LiteralPath $pidFile -Raw).Trim()
    if ($existingPidText -match "^\d+$") {
        $existingPid = [int]$existingPidText
        $existingProcess = Get-CimInstance Win32_Process -Filter "ProcessId = $existingPid" -ErrorAction SilentlyContinue
        if ($null -ne $existingProcess -and $existingProcess.CommandLine -like "*run_sdss_dr16q_v4_cross_catalogue_pipeline_2026-07-18.py*") {
            Write-Host "The SDSS DR16Q v4 pipeline is already running."
            Write-Host "PID: $existingPid"
            Write-Host "stdout log: $stdoutLog"
            Write-Host "stderr log: $stderrLog"
            exit 0
        }
    }
    Remove-Item -LiteralPath $pidFile -Force
}

$python = (Get-Command python.exe -ErrorAction Stop).Source
$process = Start-Process `
    -FilePath $python `
    -ArgumentList @("-u", $pipeline) `
    -WorkingDirectory $root `
    -WindowStyle Hidden `
    -RedirectStandardOutput $stdoutLog `
    -RedirectStandardError $stderrLog `
    -PassThru

Set-Content -LiteralPath $pidFile -Value $process.Id -Encoding ascii

Write-Host "Started preregistered SDSS DR16Q v4 cross-catalogue pipeline."
Write-Host "PID: $($process.Id)"
Write-Host "stdout log: $stdoutLog"
Write-Host "stderr log: $stderrLog"
Write-Host ""
Write-Host "Tail progress with:"
Write-Host "  Get-Content -Wait -Tail 60 `"$stdoutLog`""
