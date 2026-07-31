[CmdletBinding()]
param()

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$rustRoot = Join-Path $repoRoot 'bindings\rust'
$sourceRoot = Join-Path $rustRoot 'src'
$manifestPath = Join-Path $rustRoot 'Cargo.toml'
$buildScriptPath = Join-Path $rustRoot 'build.rs'
$errors = [System.Collections.Generic.List[string]]::new()

$forbiddenSourcePatterns = @(
    'threadpool',
    'ThreadPool',
    'std::thread',
    'mpsc',
    'sync_channel',
    'Barrier',
    'available_parallelism',
    'num_cpus',
    'joined_execute',
    'da_pool'
)

foreach ($file in @(Get-ChildItem -LiteralPath $sourceRoot -Filter '*.rs' -File))
{
    foreach ($pattern in $forbiddenSourcePatterns)
    {
        foreach ($match in @(Select-String -LiteralPath $file.FullName -Pattern $pattern -SimpleMatch))
        {
            $isCallerThreadTest =
                $file.Name -eq 'pairing.rs' -and
                $match.Line.Trim() -eq 'let pairing = std::thread::spawn(move || {'
            if (-not $isCallerThreadTest)
            {
                $errors.Add("$($file.Name):$($match.LineNumber) contains forbidden internal-threading token '$pattern'.")
            }
        }
    }
}

$manifestText = Get-Content -LiteralPath $manifestPath -Raw
foreach ($dependency in @('threadpool', 'rayon', 'crossbeam', 'num_cpus'))
{
    if ($manifestText -match "(?m)^\s*$dependency\s*=")
    {
        $errors.Add("Cargo.toml contains scheduler dependency '$dependency'.")
    }
}
if ($manifestText -match '(?m)^\s*no-threads\s*=')
{
    $errors.Add("Cargo.toml still exposes the obsolete 'no-threads' feature.")
}

$buildScriptText = Get-Content -LiteralPath $buildScriptPath -Raw
if (-not $buildScriptText.Contains('"-DBLST_PARALLEL_MSM=OFF"'))
{
    $errors.Add('build.rs does not explicitly disable native BLST parallel MSM.')
}

if ($errors.Count -gt 0)
{
    Write-Error ("Rust threading policy failed with {0} error(s):`n - {1}" -f
        $errors.Count, ($errors -join "`n - "))
    exit 1
}

Write-Output 'Rust threading policy passed: no internal scheduler or production thread creation found.'
