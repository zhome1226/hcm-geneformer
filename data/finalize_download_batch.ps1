param(
    [string]$LocalDirectory = 'D:\VirtualCell\sources\expansion',
    [string]$Repository = 'D:\VirtualCell\xiaoe-internship'
)
# Finish this one-off batch only. Publish text reports if the checkout remains
# exactly at the launch commit and clean; otherwise leave reports for review.
$ErrorActionPreference = 'Stop'
$completion = Join-Path $LocalDirectory 'final\COMPLETION.md'
$expectedFile = Join-Path $LocalDirectory 'expected_git_head.txt'
while (-not (Test-Path -LiteralPath $completion)) { Start-Sleep -Seconds 60 }
$sync = @{status='reports_saved';updated=(Get-Date -Format o)}
try {
    if (-not (Test-Path -LiteralPath $expectedFile)) { throw 'Missing expected launch commit' }
    $expected = (Get-Content -LiteralPath $expectedFile -Raw).Trim()
    $head = (& git -C $Repository rev-parse HEAD).Trim()
    $branch = (& git -C $Repository branch --show-current).Trim()
    $origin = (& git -C $Repository remote get-url origin).Trim()
    $dirty = & git -C $Repository status --porcelain
    if ($head -ne $expected -or $branch -ne 'virtual-cell' -or $dirty) {
        throw 'Checkout changed since launch; automatic publishing skipped'
    }
    if ($origin -ne 'https://github.com/zhome1226/xiaoe-internship.git') {
        throw 'Origin changed; automatic publishing skipped'
    }
    $paths = @('data/expansion/MANIFEST.tsv','data/expansion/MD5SUMS',
               'data/expansion/STATUS.tsv','data/expansion/STATUS.json')
    foreach ($path in $paths) {
        Copy-Item -LiteralPath (Join-Path $LocalDirectory ('final\'+(Split-Path $path -Leaf))) `
            -Destination (Join-Path $Repository $path)
    }
    $state = Get-Content -LiteralPath (Join-Path $LocalDirectory 'final\STATUS.json') -Raw | ConvertFrom-Json
    $message = if ($state.status -eq 'complete') {
        'data: record verified public dataset download batch'
    } else {
        'data: record partial public dataset downloads'
    }
    & git -C $Repository add -- @paths
    if ($LASTEXITCODE -ne 0) { throw 'Git add failed' }
    & git -C $Repository diff --cached --check
    if ($LASTEXITCODE -ne 0) { throw 'Text checks failed' }
    & git -C $Repository -c 'user.name=Zhang Shuo' `
        -c 'user.email=111625473+zhome1226@users.noreply.github.com' commit --only -m $message -- @paths
    if ($LASTEXITCODE -ne 0) { throw 'Git commit failed' }
    $sync.commit = (& git -C $Repository rev-parse HEAD).Trim()
    $env:GIT_TERMINAL_PROMPT = '0'
    & git -C $Repository push origin HEAD:virtual-cell
    if ($LASTEXITCODE -ne 0) { throw 'Git push failed; commit retained locally' }
    $remote = & git -C $Repository ls-remote origin refs/heads/virtual-cell
    if ($LASTEXITCODE -ne 0 -or ($remote -split '\s+')[0] -ne $sync.commit) {
        throw 'Remote commit confirmation failed'
    }
    $sync.status = 'pushed'
    $sync.download_status = $state.status
} catch {
    $sync.status = 'needs_review'
    $sync.reason = $_.Exception.Message
}
$sync.updated = Get-Date -Format o
$sync | ConvertTo-Json -Depth 6 | Set-Content -Encoding utf8 `
    (Join-Path $LocalDirectory 'final\SYNC_STATUS.json')
