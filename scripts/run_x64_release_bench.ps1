$ErrorActionPreference = 'Stop'

$vsCandidates = @(
    $env:VSINSTALLDIR,
    'C:\Program Files\Microsoft Visual Studio\2022\Community',
    'C:\Program Files\Microsoft Visual Studio\18\Insiders'
) | Where-Object { $_ -and (Test-Path $_) }

if (-not $vsCandidates) {
    throw 'Unable to locate a Visual Studio installation (checked VSINSTALLDIR, VS 2022 Community, and VS 18 Insiders).'
}

$vsInstall = $vsCandidates[0].TrimEnd('\')
Write-Host "Using Visual Studio install: $vsInstall"
$cmake     = Join-Path $vsInstall 'Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe'
$buildDir  = 'out\build\x64-Release'
$benchExe  = Join-Path $buildDir 'benchmark\blst_benchmarks.exe'
$logDir    = 'out\logs'
$logFile   = Join-Path $logDir 'x64_release_bench.txt'
$jsonFile  = Join-Path $logDir 'x64_release_bench.json'
$cudaRoot  = 'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.0'
$defaultBenchFilter = '(bench_(g1_fft_cpu/256/.*|gt_ifft_cpu/256/.*|g1_fft_gpu/256/.*|gt_ifft_gpu/256/.*|middle_product_basic/128/.*|middle_product_fft/128/.*))'
$benchFilter = if ($env:BLST_BENCH_FILTER) { $env:BLST_BENCH_FILTER } else { $defaultBenchFilter }
$benchMinTime = if ($env:BLST_BENCH_MIN_TIME) { $env:BLST_BENCH_MIN_TIME } else { '20ms' }

if (-not (Test-Path $logDir)) {
    New-Item -ItemType Directory -Path $logDir | Out-Null
}

$env:PATH = @(
    (Join-Path $vsInstall 'Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja')
    (Join-Path $vsInstall 'Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin')
    $env:PATH
) -join ';'

$msvcToolsRoot = Join-Path $vsInstall 'VC\Tools\MSVC'
$msvcVersionDir = Get-ChildItem -Directory $msvcToolsRoot | Sort-Object Name -Descending | Select-Object -First 1
if (-not $msvcVersionDir) {
    throw "No MSVC toolset directory found under $msvcToolsRoot"
}
$toolBin = Join-Path $msvcVersionDir.FullName 'bin\Hostx64\x64'
$clExe   = Join-Path $toolBin 'cl.exe'
$mlExe   = Join-Path $toolBin 'ml64.exe'
$linkExe = Join-Path $toolBin 'link.exe'
$libExe  = Join-Path $toolBin 'lib.exe'
foreach ($tool in @($clExe, $mlExe, $linkExe, $libExe)) {
    if (-not (Test-Path $tool)) {
        throw "Required MSVC tool not found at $tool"
    }
}
$compilerArgs = @(
    "-D CMAKE_C_COMPILER=`"$clExe`"",
    "-D CMAKE_CXX_COMPILER=`"$clExe`"",
    "-D CMAKE_ASM_MASM_COMPILER=`"$mlExe`"",
    "-D CMAKE_LINKER=`"$linkExe`"",
    "-D CMAKE_AR=`"$libExe`""
) -join ' '

$vcInclude = Join-Path $msvcVersionDir.FullName 'include'
$vcLib     = Join-Path $msvcVersionDir.FullName 'lib\x64'
$sdkRoot   = Join-Path ${env:ProgramFiles(x86)} 'Windows Kits\10'
$sdkIncludeRoot = Join-Path $sdkRoot 'Include'
$sdkLibRoot     = Join-Path $sdkRoot 'Lib'
$sdkVersionDir = $null
if (Test-Path $sdkIncludeRoot) {
    $sdkVersionDir = Get-ChildItem -Directory $sdkIncludeRoot | Sort-Object Name -Descending | Select-Object -First 1
}
if (-not $sdkVersionDir) {
    throw "Unable to locate Windows 10 SDK include directory under $sdkIncludeRoot"
}
$sdkVersion = $sdkVersionDir.Name
$sdkIncludePaths = @(
    $vcInclude
    (Join-Path $sdkIncludeRoot "$sdkVersion\ucrt")
    (Join-Path $sdkIncludeRoot "$sdkVersion\shared")
    (Join-Path $sdkIncludeRoot "$sdkVersion\um")
    (Join-Path $sdkIncludeRoot "$sdkVersion\winrt")
    (Join-Path $sdkIncludeRoot "$sdkVersion\cppwinrt")
) | Where-Object { $_ -and (Test-Path $_) }
$sdkLibPaths = @(
    $vcLib
    (Join-Path $sdkLibRoot "$sdkVersion\ucrt\x64")
    (Join-Path $sdkLibRoot "$sdkVersion\um\x64")
) | Where-Object { $_ -and (Test-Path $_) }

$existingInclude = if ($env:INCLUDE) { $env:INCLUDE } else { '' }
$existingLib = if ($env:LIB) { $env:LIB } else { '' }
$env:INCLUDE = (@($sdkIncludePaths) + @($existingInclude)) -join ';'
$env:LIB     = (@($sdkLibPaths) + @($existingLib)) -join ';'

$vcvars = Join-Path $vsInstall 'VC\Auxiliary\Build\vcvars64.bat'

function Invoke-Cmd {
    param([string]$Command, [string]$ErrorMessage)
    $psi = New-Object System.Diagnostics.ProcessStartInfo -Property @{
        FileName        = 'cmd.exe'
        Arguments       = "/c `"$Command`""
        UseShellExecute = $false
    }
    $process = [System.Diagnostics.Process]::Start($psi)
    $process.WaitForExit()
    if ($process.ExitCode -ne 0) {
        throw "$ErrorMessage (exit code $($process.ExitCode))"
    }
}

Write-Host '[1/3] Configure preset x64-Release'
Invoke-Cmd "call `"$vcvars`" x64 && `"$cmake`" --preset x64-Release $compilerArgs -D CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES=`"$($cudaRoot)\include;$($cudaRoot)\include\cccl`"" 'cmake configure failed'

Write-Host '[2/3] Build blst_benchmarks (RelWithDebInfo)'
Invoke-Cmd "call `"$vcvars`" x64 && `"$cmake`" --build `"$buildDir`" --target blst_benchmarks --config RelWithDebInfo" 'cmake build failed'

if (-not (Test-Path $benchExe)) {
    throw "Benchmark executable not found at $benchExe"
}

Write-Host '[3/3] Run FFT + middle-product benchmarks'
$tempSuffix = [System.IO.Path]::GetRandomFileName()
$tempLog  = Join-Path $logDir "x64_release_bench_$tempSuffix.txt"
$tempJson = Join-Path $logDir "x64_release_bench_$tempSuffix.json"
$benchCmd = "call `"$vcvars`" x64 && `"$benchExe`" --benchmark_filter=`"$benchFilter`" --benchmark_min_time=$benchMinTime --benchmark_out=`"$tempJson`" --benchmark_out_format=json"
Invoke-Cmd "$benchCmd > `"$tempLog`"" 'benchmark run failed'

if (Test-Path $tempJson) {
    Remove-Item -Force -ErrorAction SilentlyContinue $jsonFile
    Move-Item -Force -Path $tempJson -Destination $jsonFile
}
if (Test-Path $tempLog) {
    Remove-Item -Force -ErrorAction SilentlyContinue $logFile
    Move-Item -Force -Path $tempLog -Destination $logFile
}

Write-Host "Console log: $logFile"
Write-Host "JSON output: $jsonFile"
