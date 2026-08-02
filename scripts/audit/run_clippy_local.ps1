[CmdletBinding()]
param(
    [ValidateSet('default', 'full', 'portable', 'avx', 'full-serde', 'all')]
    [string]$Mode = 'all'
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$rustRoot = Join-Path $repoRoot 'bindings\rust'
$commit = (git -C $repoRoot rev-parse HEAD).Trim()
$shortCommit = $commit.Substring(0, 12)
$timestamp = (Get-Date).ToUniversalTime().ToString('yyyyMMddTHHmmssZ')
$runId = "$shortCommit-$timestamp"
$resultsRoot = Join-Path $repoRoot "out\audit\clippy\$runId"

if (Test-Path -LiteralPath $resultsRoot)
{
    throw "Refusing to overwrite an existing Clippy run: $runId"
}

$profiles = [ordered]@{
    'default' = ''
    'full' = 'full-rust-bindings'
    'portable' = 'portable'
    'avx' = 'avx'
    'full-serde' = 'full-rust-bindings,serde-secret'
}

if ($Mode -eq 'all')
{
    $selectedProfiles = @($profiles.Keys)
}
else
{
    $selectedProfiles = @($Mode)
}

New-Item -ItemType Directory -Path $resultsRoot | Out-Null
git -C $repoRoot status --porcelain=v1 --untracked-files=all |
    Set-Content -LiteralPath (Join-Path $resultsRoot 'working-tree-status.txt')

$rustcVersion = (& rustc --version --verbose) -join "`n"
$cargoVersion = (& cargo --version) -join "`n"
$clippyVersion = (& cargo clippy --version) -join "`n"
@(
    "commit=$commit"
    "timestamp_utc=$timestamp"
    "mode=$Mode"
    'cargo_profile=dev'
    'targets=library and build script'
    'lints=clippy::all,clippy::pedantic,clippy::nursery'
    "rustc_version=$($rustcVersion -replace "`r?`n", '; ')"
    "cargo_version=$cargoVersion"
    "clippy_version=$clippyVersion"
) | Set-Content -LiteralPath (Join-Path $resultsRoot 'run-metadata.txt')

$env:CARGO_TERM_COLOR = 'never'
$summary = [System.Collections.Generic.List[object]]::new()

Push-Location $rustRoot
try
{
    foreach ($profile in $selectedProfiles)
    {
        $features = $profiles[$profile]
        $profileRoot = Join-Path $resultsRoot $profile
        $jsonLog = Join-Path $profileRoot 'clippy.jsonl'
        $stderrLog = Join-Path $profileRoot 'cargo.stderr.log'
        New-Item -ItemType Directory -Path $profileRoot | Out-Null

        $cargoArgs = @(
            'clippy'
            '--lib'
            '--no-default-features'
            '--message-format=json'
        )
        if ($features.Length -ne 0)
        {
            $cargoArgs += @('--features', $features)
        }
        $cargoArgs += @(
            '--'
            '-W'
            'clippy::all'
            '-W'
            'clippy::pedantic'
            '-W'
            'clippy::nursery'
        )

        Write-Output "Running Clippy profile '$profile' (features='$features')"
        $savedErrorActionPreference = $ErrorActionPreference
        $ErrorActionPreference = 'Continue'
        & cargo @cargoArgs 2> $stderrLog |
            Tee-Object -FilePath $jsonLog |
            Out-Null
        $cargoExitCode = $LASTEXITCODE
        $ErrorActionPreference = $savedErrorActionPreference
        if ($cargoExitCode -ne 0)
        {
            throw "Clippy profile '$profile' failed with exit code $cargoExitCode"
        }

        $diagnostics = foreach ($line in Get-Content -LiteralPath $jsonLog)
        {
            try
            {
                $item = $line | ConvertFrom-Json
            }
            catch
            {
                continue
            }
            if ($item.reason -ne 'compiler-message' -or
                $item.message.level -ne 'warning')
            {
                continue
            }
            $span = $item.message.spans |
                Where-Object is_primary |
                Select-Object -First 1
            if ($null -eq $span)
            {
                continue
            }
            [pscustomobject]@{
                Key = "$($item.message.code.code)|$($span.file_name)|$($span.line_start)|$($span.column_start)"
                Lint = $item.message.code.code
                File = $span.file_name
                Line = $span.line_start
                Column = $span.column_start
                Message = $item.message.message
            }
        }

        $uniqueDiagnostics = @($diagnostics | Sort-Object Key -Unique)
        $uniqueDiagnostics |
            Select-Object Lint, File, Line, Column, Message |
            Export-Csv -NoTypeInformation -LiteralPath (Join-Path $profileRoot 'diagnostics.csv')
        $lintCounts = @($uniqueDiagnostics | Group-Object Lint | Sort-Object Count -Descending)
        $lintCounts |
            ForEach-Object {
                [pscustomobject]@{ Count = $_.Count; Lint = $_.Name }
            } |
            Export-Csv -NoTypeInformation -LiteralPath (Join-Path $profileRoot 'lint-counts.csv')

        $summary.Add([pscustomobject]@{
            Profile = $profile
            Features = $features
            RawWarnings = @($diagnostics).Count
            CanonicalWarnings = $uniqueDiagnostics.Count
        })
        Write-Output "  canonical warnings: $($uniqueDiagnostics.Count)"
    }
}
finally
{
    Pop-Location
}

$summary | Export-Csv -NoTypeInformation -LiteralPath (Join-Path $resultsRoot 'summary.csv')
Write-Output "Clippy results: $resultsRoot"
