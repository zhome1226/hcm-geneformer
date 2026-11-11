param(
    [string]$LocalDirectory = 'D:\VirtualCell\sources\expansion',
    [string]$Remote = 'zhome@10.158.203.133',
    [int]$ExistingTunnelProcessId = 0
)
# One-off transport helper for the currently running download batch.
# Keeps SSH forwarding alive; on completion copies factual reports, then exits.
$ErrorActionPreference = 'Stop'
New-Item -ItemType Directory -Force -Path $LocalDirectory | Out-Null
$ssh = (Get-Command ssh.exe).Source
$scp = (Get-Command scp.exe).Source
$common = @('-o','BatchMode=yes','-o','ConnectTimeout=15','-o','IPQoS=none',
            '-o','ServerAliveInterval=10','-o','ServerAliveCountMax=3')
$stopFile = Join-Path $LocalDirectory 'STOP_DOWNLOAD_MONITOR'
$tunnel = $null
if ($ExistingTunnelProcessId -gt 0) {
    try { $tunnel = Get-Process -Id $ExistingTunnelProcessId } catch { $tunnel = $null }
}
try {
    while (-not (Test-Path -LiteralPath $stopFile)) {
        if ($null -eq $tunnel -or $tunnel.HasExited) {
            $arguments = @('-N') + $common + @('-o','ExitOnForwardFailure=yes',
                '-R','127.0.0.1:39843:127.0.0.1:7890',$Remote)
            $tunnel = Start-Process -FilePath $ssh -ArgumentList $arguments -PassThru `
                -WindowStyle Hidden -RedirectStandardError (Join-Path $LocalDirectory 'tunnel.stderr.log')
            Start-Sleep -Seconds 3
            if ($tunnel.HasExited) { Start-Sleep -Seconds 10; continue }
        }
        try {
            $queryArguments = @('-n','-T') + $common + @('-o','ClearAllForwardings=yes',
                $Remote,'cat ~/xiaoe-internship-data/_download_expansion/STATUS.json')
            $rawPath = Join-Path $LocalDirectory 'status.raw.json'
            $query = Start-Process -FilePath $ssh -ArgumentList $queryArguments -PassThru `
                -WindowStyle Hidden -RedirectStandardOutput $rawPath `
                -RedirectStandardError (Join-Path $LocalDirectory 'status.stderr.log')
            if (-not $query.WaitForExit(30000)) {
                Stop-Process -Id $query.Id
                throw 'Status SSH timed out after 30 seconds'
            }
            if ($query.ExitCode -ne 0) { throw 'Status SSH failed' }
            $raw = Get-Content -LiteralPath $rawPath -Raw
            $state = $raw | ConvertFrom-Json
            $raw | Set-Content -Encoding utf8 (Join-Path $LocalDirectory 'LIVE_STATUS.json')
            if ($state.status -in @('complete','partial_failure')) {
                $destination = Join-Path $LocalDirectory 'final'
                New-Item -ItemType Directory -Force -Path $destination | Out-Null
                $remoteBase = $Remote + ':/home/zhome/xiaoe-internship-data/_download_expansion/'
                & $scp @common ($remoteBase+'MANIFEST.tsv') ($remoteBase+'MD5SUMS') `
                    ($remoteBase+'STATUS.tsv') ($remoteBase+'STATUS.json') $destination
                if ($LASTEXITCODE -ne 0) { throw 'Final report copy failed' }
                @(
                    '# Download batch finished'
                    ('Status: ' + $state.status)
                    ('Verified files: ' + $state.files_verified + '/' + $state.files_planned)
                    ('Planned bytes: ' + $state.planned_bytes)
                    ('UTC: ' + $state.updated_utc)
                    'See STATUS.json for errors; MANIFEST.tsv contains verified files only.'
                    'No expression analysis, training, or Git publishing was performed by this helper.'
                ) | Set-Content -Encoding utf8 (Join-Path $destination 'COMPLETION.md')
                break
            }
        } catch {
            (Get-Date -Format o) + ' ' + $_.Exception.Message | Add-Content `
                (Join-Path $LocalDirectory 'monitor.stderr.log')
        }
        Start-Sleep -Seconds 60
    }
} finally {
    if ($null -ne $tunnel -and -not $tunnel.HasExited) {
        Stop-Process -Id $tunnel.Id
    }
}
