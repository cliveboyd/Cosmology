$ErrorActionPreference = "Stop"

$TaskRoot = "C:\Users\clive\Documents\Codex\2026-07-14\ok-continu"
$Pipeline = Join-Path $TaskRoot "outputs\run_su2_quaia_targeted_stability_mocks.ps1"
$LogOut = Join-Path $TaskRoot "outputs\su2_quaia_targeted_stability_mocks_stdout.log"
$LogErr = Join-Path $TaskRoot "outputs\su2_quaia_targeted_stability_mocks_stderr.log"

if (-not (Test-Path -LiteralPath $Pipeline)) {
  throw "Missing pipeline script: $Pipeline"
}

if (Test-Path -LiteralPath $LogOut) {
  $stamp = Get-Date -Format "yyyyMMdd_HHmmss"
  Rename-Item -LiteralPath $LogOut -NewName "su2_quaia_targeted_stability_mocks_stdout_$stamp.log"
}
if (Test-Path -LiteralPath $LogErr) {
  $stamp = Get-Date -Format "yyyyMMdd_HHmmss"
  Rename-Item -LiteralPath $LogErr -NewName "su2_quaia_targeted_stability_mocks_stderr_$stamp.log"
}

$Args = @(
  "-NoProfile",
  "-ExecutionPolicy", "Bypass",
  "-File", $Pipeline
)

$Process = Start-Process `
  -FilePath "powershell.exe" `
  -ArgumentList $Args `
  -WorkingDirectory $TaskRoot `
  -WindowStyle Hidden `
  -RedirectStandardOutput $LogOut `
  -RedirectStandardError $LogErr `
  -PassThru

Write-Host "Started SU2/Quaia targeted stability mock percentile audit."
Write-Host "PID: $($Process.Id)"
Write-Host "stdout log: $LogOut"
Write-Host "stderr log: $LogErr"
Write-Host ""
Write-Host "Tail the log with:"
Write-Host "  Get-Content -Wait -Tail 80 `"$LogOut`""
