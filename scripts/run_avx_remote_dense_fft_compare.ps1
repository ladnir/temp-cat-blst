param(
    [string]$RemoteBase = "/tmp/avx-dense-fft",
    [string]$ClangRoot = "/opt/llvm-21/bin",
    [int]$PinCore = 0,
    [int]$Repetitions = 5,
    [string]$MinTime = "0.05s",
    [string]$SizesRegex = "2|4|8|16|32"
)

$ErrorActionPreference = "Stop"

$rootDir = Split-Path -Parent $PSScriptRoot
$connectScript = Join-Path $env:USERPROFILE ".codex\skills\connect-to-aws\scripts\connect-ec2.ps1"
$plink = "C:\Users\peter\OneDrive\tools\plink.exe"
$key = "C:\Users\peter\OneDrive\keys\aws_2_19_26.ppk"
$remoteHost = "13.58.237.236"
$user = "ec2-user"
$hostKey = "ssh-ed25519 255 79:1d:89:3c:8a:26:00:c7:6a:fb:b2:1c:24:e3:c8:6e"
$runId = [DateTime]::UtcNow.ToString("yyyyMMddTHHmmssZ")
$remoteRunDir = "$RemoteBase/$runId"
$localReportDir = Join-Path $rootDir "out\avx-dense-fft\$runId"

New-Item -ItemType Directory -Force -Path $localReportDir | Out-Null

function Invoke-Remote {
    param([string]$Command)

    powershell -NoProfile -NonInteractive -ExecutionPolicy Bypass -File $connectScript -Command $Command
}

function Upload-File {
    param(
        [string]$Source,
        [string]$Destination
    )

    $cmd = "`"$plink`" -batch -ssh -l $user -i `"$key`" -hostkey `"$hostKey`" $remoteHost `"cat > '$Destination'`" < `"$Source`""
    cmd.exe /c $cmd | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw "upload failed for $Source"
    }
}

$tmpDir = Join-Path ([System.IO.Path]::GetTempPath()) ("blst-dense-fft-" + [System.Guid]::NewGuid().ToString("N"))
$snapshotDir = Join-Path $tmpDir "snapshot"
$snapshotTar = Join-Path $tmpDir "blst_current.tar"
$snapshotArchive = Join-Path $tmpDir "blst_current.tar.gz"
New-Item -ItemType Directory -Force -Path $snapshotDir | Out-Null

try {
    Push-Location $rootDir

    & git archive --format=tar -o $snapshotTar HEAD
    tar -xf $snapshotTar -C $snapshotDir

    $changedPaths = & git diff --name-only HEAD
    foreach ($path in $changedPaths) {
        if ([string]::IsNullOrWhiteSpace($path)) {
            continue
        }

        $sourcePath = Join-Path $rootDir $path
        $destPath = Join-Path $snapshotDir $path
        if (Test-Path $sourcePath -PathType Leaf) {
            New-Item -ItemType Directory -Force -Path (Split-Path -Parent $destPath) | Out-Null
            Copy-Item -LiteralPath $sourcePath -Destination $destPath -Force
        } elseif (Test-Path $destPath) {
            Remove-Item -LiteralPath $destPath -Force
        }
    }

    $untrackedPaths = & git ls-files --others --exclude-standard
    foreach ($path in $untrackedPaths) {
        if ([string]::IsNullOrWhiteSpace($path)) {
            continue
        }

        $sourcePath = Join-Path $rootDir $path
        if (-not (Test-Path $sourcePath -PathType Leaf)) {
            continue
        }

        $destPath = Join-Path $snapshotDir $path
        New-Item -ItemType Directory -Force -Path (Split-Path -Parent $destPath) | Out-Null
        Copy-Item -LiteralPath $sourcePath -Destination $destPath -Force
    }

    tar -czf $snapshotArchive -C $snapshotDir .
}
finally {
    Pop-Location
}

Write-Host "remote run dir: $remoteRunDir"
Write-Host "local report dir: $localReportDir"

Invoke-Remote "mkdir -p '$remoteRunDir/archives'"
Upload-File $snapshotArchive "$remoteRunDir/archives/blst_current.tar.gz"
Invoke-Remote "rm -rf '$remoteRunDir/blst_current' && mkdir -p '$remoteRunDir/blst_current' && tar -xzf '$remoteRunDir/archives/blst_current.tar.gz' -C '$remoteRunDir/blst_current'"

$configureAndBuild = @"
cd '$remoteRunDir/blst_current' && \
cmake -S . -B out/dense-fft-clang21 -G Ninja \
  -DCMAKE_C_COMPILER='$ClangRoot/clang' \
  -DCMAKE_CXX_COMPILER='$ClangRoot/clang++' \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_FLAGS_RELEASE='-O3 -march=native -DNDEBUG' \
  -DCMAKE_CXX_FLAGS_RELEASE='-O3 -march=native -DNDEBUG' \
  -DBLST_PORTABLE=OFF -DBLST_ASM=ON -DBLST_ADX=ON \
  -DBLST_AVX512=ON -DBLST_DEFAULT_AVX=ON -DBLST_AVX512_PORTABLE=OFF \
  -DBLST_BENCHMARK=ON -DBLST_TEST=ON -DBLST_CUDA=OFF >/dev/null && \
cmake --build out/dense-fft-clang21 --target blst_benchmarks blst_tests -j2 >/dev/null
"@
Invoke-Remote $configureAndBuild

Invoke-Remote "cd '$remoteRunDir/blst_current' && ./out/dense-fft-clang21/test/blst_tests --gtest_filter='Fft.AvxDenseTransformMatchesFft' >/dev/null"

$fftJsonPath = Join-Path $localReportDir "avx_fft.json"
$denseJsonPath = Join-Path $localReportDir "avx_dense.json"

Invoke-Remote "cd '$remoteRunDir/blst_current/out/dense-fft-clang21/benchmark' && taskset -c $PinCore ./blst_benchmarks --benchmark_filter='^bench_avx_g1_fft_cpu/($SizesRegex)/real_time$' --benchmark_repetitions=$Repetitions --benchmark_report_aggregates_only=true --benchmark_min_time=$MinTime --benchmark_format=json" | Out-File -Encoding ascii $fftJsonPath
Invoke-Remote "cd '$remoteRunDir/blst_current/out/dense-fft-clang21/benchmark' && taskset -c $PinCore ./blst_benchmarks --benchmark_filter='^bench_avx_g1_dense_transform_cpu/($SizesRegex)/real_time$' --benchmark_repetitions=$Repetitions --benchmark_report_aggregates_only=true --benchmark_min_time=$MinTime --benchmark_format=json" | Out-File -Encoding ascii $denseJsonPath

$summaryPath = Join-Path $localReportDir "summary.txt"

function Load-BenchmarkMap {
    param([string]$Path)

    $text = Get-Content -Path $Path -Raw
    $jsonStart = $text.IndexOf("{")
    if ($jsonStart -lt 0) {
        throw "benchmark json not found in $Path"
    }

    $payload = $text.Substring($jsonStart) | ConvertFrom-Json
    $map = @{}

    $withMedian = @($payload.benchmarks | Where-Object { $_.aggregate_name -eq "median" })
    if ($withMedian.Count -gt 0) {
        foreach ($bench in $withMedian) {
            $map[$bench.name] = @{
                time = [double]$bench.cpu_time
                unit = [string]$bench.time_unit
            }
        }
        return $map
    }

    foreach ($bench in $payload.benchmarks) {
        $map[$bench.name] = @{
            time = [double]$bench.cpu_time
            unit = [string]$bench.time_unit
        }
    }
    return $map
}

$fftMap = Load-BenchmarkMap $fftJsonPath
$denseMap = Load-BenchmarkMap $denseJsonPath
$lines = New-Object System.Collections.Generic.List[string]

foreach ($fftName in ($fftMap.Keys | Sort-Object)) {
    if ($fftName -notmatch "/(\d+)(/|$)") {
        continue
    }

    $size = $Matches[1]
    $denseName = "bench_avx_g1_dense_transform_cpu/$size"
    if (-not $denseMap.ContainsKey($denseName)) {
        continue
    }

    $fftTime = $fftMap[$fftName].time
    $denseTime = $denseMap[$denseName].time
    $unit = $fftMap[$fftName].unit
    $ratio = $denseTime / $fftTime
    $drift = (($denseTime - $fftTime) / $fftTime) * 100.0
    $lines.Add(("N={0}: fft {1:N3} {2}, dense {3:N3} {2}, drift {4:+0.000;-0.000}%, dense/fft {5:N3}x" -f $size, $fftTime, $unit, $denseTime, $drift, $ratio))
}

$lines | Set-Content -Encoding ascii $summaryPath
$lines | ForEach-Object { Write-Host $_ }

Write-Host "summary written to $summaryPath"

Remove-Item -LiteralPath $tmpDir -Recurse -Force
