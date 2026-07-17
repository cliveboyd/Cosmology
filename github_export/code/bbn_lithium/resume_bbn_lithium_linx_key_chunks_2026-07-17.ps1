param(
    [string]$Python = "C:\Users\clive\Documents\Cosmology\plamb_runs\envs\linx_py311\python.exe",
    [string]$Programme = "C:\Users\clive\Documents\Cosmology\github_export\code\bbn_lithium\analyze_bbn_lithium_linx_fr_network_2026-07-16.py",
    [string]$OutDir = "C:\Users\clive\Documents\Cosmology\plamb_runs\diagnostics\bbn_lithium_linx_key_fr_network_overnight_20260716",
    [string]$LogDir = "C:\Users\clive\Documents\Codex\2026-07-14\ok-continu\outputs"
)

$ErrorActionPreference = "Stop"
$pointsPath = Join-Path $OutDir "bbn_lithium_linx_fr_network_points.csv"
$stages = @(23000, 24000, 25000, 26000, 27000, 28000, 29000, 30000, 30503)

function Get-UniqueStatus {
    $rows = Import-Csv -LiteralPath $pointsPath
    $okIds = [Collections.Generic.HashSet[string]]::new([StringComparer]::Ordinal)
    $failedIds = [Collections.Generic.HashSet[string]]::new([StringComparer]::Ordinal)
    foreach ($row in $rows) {
        if ($row.status -eq "ok") {
            [void]$okIds.Add($row.point_id)
        }
    }
    foreach ($row in $rows) {
        if ($row.status -eq "failed" -and -not $okIds.Contains($row.point_id)) {
            [void]$failedIds.Add($row.point_id)
        }
    }
    [pscustomobject]@{
        Ok         = $okIds.Count
        Unresolved = $failedIds.Count
    }
}

New-Item -ItemType Directory -Force -Path $LogDir | Out-Null
Write-Output "[$(Get-Date -Format o)] Starting chunked LINX resume"
Write-Output "Python: $Python"
Write-Output "Programme: $Programme"
Write-Output "Output: $OutDir"
Write-Output "Initial status: $((Get-UniqueStatus) | ConvertTo-Json -Compress)"

foreach ($stage in $stages) {
    $stdout = Join-Path $LogDir ("bbn_lithium_linx_resume_{0}_stdout.log" -f $stage)
    $stderr = Join-Path $LogDir ("bbn_lithium_linx_resume_{0}_stderr.log" -f $stage)
    Write-Output "[$(Get-Date -Format o)] Stage max_points=$stage"

    & $Python $Programme `
        --network key `
        --preset overnight `
        --resume `
        --outdir $OutDir `
        --max-points $stage `
        --n-random 30000 `
        --seed 20260716 `
        --rtol 1e-5 `
        --atol 1e-9 `
        --sampling-key 150 `
        --summary-every 250 `
        1>> $stdout 2>> $stderr

    if ($LASTEXITCODE -ne 0) {
        throw "LINX stage $stage exited with code $LASTEXITCODE. See $stderr"
    }
    $status = Get-UniqueStatus
    Write-Output "[$(Get-Date -Format o)] Stage $stage finished: ok=$($status.Ok) unresolved=$($status.Unresolved)"
}

# One fresh final retry resolves any point that failed late in the last chunk.
$finalStatus = Get-UniqueStatus
if ($finalStatus.Ok -lt 30503) {
    $stdout = Join-Path $LogDir "bbn_lithium_linx_resume_30503_retry_stdout.log"
    $stderr = Join-Path $LogDir "bbn_lithium_linx_resume_30503_retry_stderr.log"
    Write-Output "[$(Get-Date -Format o)] Final retry: ok=$($finalStatus.Ok) unresolved=$($finalStatus.Unresolved)"
    & $Python $Programme `
        --network key `
        --preset overnight `
        --resume `
        --outdir $OutDir `
        --max-points 30503 `
        --n-random 30000 `
        --seed 20260716 `
        --rtol 1e-5 `
        --atol 1e-9 `
        --sampling-key 150 `
        --summary-every 250 `
        1>> $stdout 2>> $stderr
    if ($LASTEXITCODE -ne 0) {
        throw "Final LINX retry exited with code $LASTEXITCODE. See $stderr"
    }
}

$finalStatus = Get-UniqueStatus
Write-Output "[$(Get-Date -Format o)] Final status: ok=$($finalStatus.Ok) unresolved=$($finalStatus.Unresolved)"
if ($finalStatus.Ok -ne 30503 -or $finalStatus.Unresolved -ne 0) {
    throw "LINX resume finished with unresolved points."
}
