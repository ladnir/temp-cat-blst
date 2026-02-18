param(
    [string]$OutTex = "docs/planning/fft_middle_product_report.tex"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Load-Json([string]$Path) {
    if (-not (Test-Path $Path)) {
        throw "Missing benchmark file: $Path"
    }
    return (Get-Content $Path -Raw | ConvertFrom-Json)
}

function Get-MeanRows([string]$Path) {
    $rows = @()
    $j = Load-Json $Path
    foreach ($b in $j.benchmarks) {
        if ($b.aggregate_name -eq "mean") {
            $rows += $b
        }
    }
    return $rows
}

function Get-IterationRows([string]$Path) {
    $rows = @()
    $j = Load-Json $Path
    foreach ($b in $j.benchmarks) {
        if ($b.run_type -eq "iteration") {
            $rows += $b
        }
    }
    return $rows
}

function Get-AggregateRows([string]$Path, [string]$AggregateName) {
    $rows = @()
    $j = Load-Json $Path
    foreach ($b in $j.benchmarks) {
        if ($b.aggregate_name -eq $AggregateName) {
            $rows += $b
        }
    }
    return $rows
}

function Get-MiddleRows([string]$BaseName) {
    $r5 = "docs/planning/$BaseName.r5.json"
    $base = "docs/planning/$BaseName.json"
    if (Test-Path $r5) {
        return Get-AggregateRows $r5 "median"
    }
    return Get-IterationRows $base
}

function Pct([double]$x) {
    return [Math]::Round($x, 1)
}

function Ms([double]$x) {
    return [Math]::Round($x, 1)
}

function Get-TotalMs($row) {
    if ($row.PSObject.Properties.Name -contains "online_total_ms") {
        return [double]$row.online_total_ms
    }
    return [double]$row.real_time
}

$g1FullRows = Get-MeanRows "docs/planning/bench_g1_fft_full_middle_shape_detailed.json"
$g1TruncRows = Get-MeanRows "docs/planning/bench_g1_fft_trunc_middle_shape_detailed.json"
$gtFullRows = Get-MeanRows "docs/planning/bench_gt_ifft_full_middle_shape_detailed.json"
$gtTruncRows = Get-MeanRows "docs/planning/bench_gt_ifft_trunc_middle_shape_detailed.json"

$basicNoTruncRows = Get-MiddleRows "bench_middle_product_basic_no_trunc_detailed"
$basicIfftTruncRows = Get-MiddleRows "bench_middle_product_basic_ifft_trunc_detailed"
$basicPrunedRows = Get-MiddleRows "bench_middle_product_basic_pruned_detailed"
$specRows = Get-MiddleRows "bench_middle_product_fft_detailed"
$specTruncRows = Get-MiddleRows "bench_middle_product_specialized_trunc_detailed"
$oneIfftRows = Get-MiddleRows "bench_middle_product_one_ifft_detailed"
$usingR5Median = (Test-Path "docs/planning/bench_middle_product_basic_no_trunc_detailed.r5.json")

$g1TruncByN = @{}
foreach ($r in $g1TruncRows) { $g1TruncByN[[int]$r.n] = $r }
$gtTruncByN = @{}
foreach ($r in $gtTruncRows) { $gtTruncByN[[int]$r.n] = $r }
$basicIfftTruncByN = @{}
foreach ($r in $basicIfftTruncRows) { $basicIfftTruncByN[[int]$r.n] = $r }
$basicPrunedByN = @{}
foreach ($r in $basicPrunedRows) { $basicPrunedByN[[int]$r.n] = $r }
$specByN = @{}
foreach ($r in $specRows) { $specByN[[int]$r.n] = $r }
$specTruncByN = @{}
foreach ($r in $specTruncRows) { $specTruncByN[[int]$r.n] = $r }
$oneIfftByN = @{}
foreach ($r in $oneIfftRows) { $oneIfftByN[[int]$r.n] = $r }

$nValues = @(64..512 | Where-Object { $_ % 32 -eq 0 })

$g1Savings = @()
$gtSavings = @()
$mpIfftTruncSavings = @()
$mpPrunedSavings = @()
$mpSpecTruncSavings = @()
$mpOneIfftVsSpecSavings = @()

$lines = New-Object System.Collections.Generic.List[string]
$lines.Add("\documentclass[11pt]{article}")
$lines.Add("\usepackage[margin=1in]{geometry}")
$lines.Add("\usepackage{booktabs}")
$lines.Add("\usepackage{longtable}")
$lines.Add("\usepackage{array}")
$lines.Add("\usepackage{siunitx}")
$lines.Add("\usepackage{hyperref}")
$lines.Add("\sisetup{round-mode=places,round-precision=1}")
$lines.Add("\begin{document}")
$lines.Add("")
$lines.Add("\title{FFT Truncation and Middle-Product Integration Report}")
$lines.Add("\author{blst benchmark sweep}")
$lines.Add("\date{2026-02-17}")
$lines.Add("\maketitle")
$lines.Add("")
$lines.Add("\section{Scope}")
$lines.Add("This report compares full and truncated FFT/iFFT kernels and their integration into middle-product flows.")
$lines.Add("Reader assumption: familiarity with FFT and middle-product algorithms.")
$lines.Add("")
$lines.Add("Compared middle-product variants:")
$lines.Add("\begin{enumerate}")
$lines.Add("\item basic FFT + basic iFFT")
$lines.Add("\item basic FFT + truncated iFFT")
$lines.Add("\item truncated FFT + truncated iFFT (pruned basic)")
$lines.Add("\item specialized middle product")
$lines.Add("\item specialized middle product + truncation")
$lines.Add("\item one-iFFT cyclic middle product")
$lines.Add("\end{enumerate}")
$lines.Add("")
$lines.Add("\section{Benchmark Setup}")
$lines.Add("\begin{itemize}")
$lines.Add("\item Range: \(n \in \{64,96,\ldots,512\}\) (step 32, non-power-of-two included).")
$lines.Add("\item Input shape: \(A\) has \(n\) coefficients in \(G_1\), \(B\) has \(2n\) coefficients in \(G_2\).")
$lines.Add("\item FFT length for basic path: \(m = 2^{\lceil \log_2(3n-1) \rceil}\).")
$lines.Add("\item Standalone kernel runs: 3 repetitions, mean reported.")
if ($usingR5Median) {
    $lines.Add("\item Integrated middle-product runs: 5 repetitions per \(n\), median aggregate reported.")
}
else {
    $lines.Add("\item Integrated middle-product runs: 1 repetition per \(n\) (long-running pairings).")
}
$lines.Add("\item Timings are real-time milliseconds from Google Benchmark JSON outputs in \texttt{docs/planning/*.json}.")
$lines.Add("\end{itemize}")
$lines.Add("")

$lines.Add("\section{Standalone Kernel Results}")
$lines.Add("\subsection{G1 Forward FFT: full vs truncated}")
$lines.Add("\begin{longtable}{r r r r r}")
$lines.Add("\toprule")
$lines.Add('$n$ & $m$ & full FFT (ms) & trunc FFT (ms) & savings (\%) \\')
$lines.Add("\midrule")
$lines.Add("\endhead")
foreach ($r in ($g1FullRows | Sort-Object { [int]$_.n })) {
    $n = [int]$r.n
    $m = [int]$r.m
    $full = [double]$r.real_time
    $tr = [double]$g1TruncByN[$n].real_time
    $sav = 100.0 * ($full - $tr) / $full
    $g1Savings += $sav
    $lines.Add(("{0} & {1} & {2} & {3} & {4} \\" -f $n, $m, (Ms $full), (Ms $tr), (Pct $sav)))
}
$lines.Add("\bottomrule")
$lines.Add("\end{longtable}")
$lines.Add("")

$lines.Add("\subsection{GT Inverse FFT: full vs truncated-range}")
$lines.Add("\begin{longtable}{r r r r r}")
$lines.Add("\toprule")
$lines.Add('$n$ & $m$ & full iFFT (ms) & trunc iFFT (ms) & savings (\%) \\')
$lines.Add("\midrule")
$lines.Add("\endhead")
foreach ($r in ($gtFullRows | Sort-Object { [int]$_.n })) {
    $n = [int]$r.n
    $m = [int]$r.m
    $full = [double]$r.real_time
    $tr = [double]$gtTruncByN[$n].real_time
    $sav = 100.0 * ($full - $tr) / $full
    $gtSavings += $sav
    $lines.Add(("{0} & {1} & {2} & {3} & {4} \\" -f $n, $m, (Ms $full), (Ms $tr), (Pct $sav)))
}
$lines.Add("\bottomrule")
$lines.Add("\end{longtable}")
$lines.Add("")

$lines.Add("\section{Integrated Middle-Product Results}")
$lines.Add("\begin{longtable}{r r r r r r r r}")
$lines.Add("\toprule")
$lines.Add('$n$ & $m$ & basic/base & basic+trunc iFFT & pruned basic & specialized & specialized+trunc & one-iFFT \\')
$lines.Add("\midrule")
$lines.Add("\endhead")
foreach ($r in ($basicNoTruncRows | Sort-Object { [int]$_.n })) {
    $n = [int]$r.n
    $m = [int]$r.m
    $base = Get-TotalMs $r
    $ifftT = Get-TotalMs $basicIfftTruncByN[$n]
    $pruned = Get-TotalMs $basicPrunedByN[$n]
    $spec = Get-TotalMs $specByN[$n]
    $specT = Get-TotalMs $specTruncByN[$n]
    $oneifft = Get-TotalMs $oneIfftByN[$n]
    $mpIfftTruncSavings += 100.0 * ($base - $ifftT) / $base
    $mpPrunedSavings += 100.0 * ($base - $pruned) / $base
    $mpSpecTruncSavings += 100.0 * ($spec - $specT) / $spec
    $mpOneIfftVsSpecSavings += 100.0 * ($spec - $oneifft) / $spec
    $lines.Add(("{0} & {1} & {2} & {3} & {4} & {5} & {6} & {7} \\" -f $n, $m, (Ms $base), (Ms $ifftT), (Ms $pruned), (Ms $spec), (Ms $specT), (Ms $oneifft)))
}
$lines.Add("\bottomrule")
$lines.Add("\end{longtable}")
$lines.Add("")

$g1Avg = ($g1Savings | Measure-Object -Average).Average
$gtAvg = ($gtSavings | Measure-Object -Average).Average
$mpIfftAvg = ($mpIfftTruncSavings | Measure-Object -Average).Average
$mpPrunedAvg = ($mpPrunedSavings | Measure-Object -Average).Average
$mpSpecAvg = ($mpSpecTruncSavings | Measure-Object -Average).Average
$mpOneIfftAvg = ($mpOneIfftVsSpecSavings | Measure-Object -Average).Average

$lines.Add("\section{Savings Summary}")
$lines.Add("\begin{itemize}")
$lines.Add(("\item Standalone forward FFT truncation average savings: {0}\%." -f (Pct $g1Avg)))
$lines.Add(("\item Standalone inverse FFT truncation average savings: {0}\%." -f (Pct $gtAvg)))
$lines.Add(("\item Integrated basic path (truncated iFFT vs basic baseline): {0}\% average savings." -f (Pct $mpIfftAvg)))
$lines.Add(("\item Integrated pruned basic path (truncated FFT+iFFT vs basic baseline): {0}\% average savings." -f (Pct $mpPrunedAvg)))
$lines.Add(("\item Specialized truncated vs specialized baseline: {0}\% average savings." -f (Pct $mpSpecAvg)))
$lines.Add(("\item One-iFFT cyclic vs specialized baseline: {0}\% average savings." -f (Pct $mpOneIfftAvg)))
$lines.Add("\end{itemize}")
$lines.Add("")

$lines.Add("\section{Interpretation}")
$lines.Add("\begin{itemize}")
$lines.Add("\item Forward truncation in \(G_1\) gives a modest positive average win, with clear size-dependent variability.")
$lines.Add("\item Standalone truncated iFFT (prefix form) is mixed in this sweep and is slightly negative on average by itself.")
$lines.Add("\item In the integrated basic middle-product flow, adding truncated iFFT is still a strong net positive.")
$lines.Add("\item The current pruned-basic implementation is not yet a net win over basic+truncated-iFFT and needs further kernel/layout tuning.")
$lines.Add("\item Specialized+trunc improves over specialized baseline and is the preferred specialized variant in this sweep.")
$lines.Add("\item One-iFFT cyclic directly captures the single-iFFT idea and should be compared head-to-head with specialized+trunc for default selection.")
$lines.Add("\item Size transitions where \(m\) doubles remain visible; non-power-of-two \(n\) values are still important for policy tuning.")
$lines.Add("\end{itemize}")
$lines.Add("")

$lines.Add("\section{Repro Commands}")
$lines.Add("\begin{verbatim}")
$lines.Add("out\\build\\x64-Release\\benchmark\\blst_benchmarks.exe --benchmark_filter=bench_g1_fft_full_middle_shape_detailed")
$lines.Add("out\\build\\x64-Release\\benchmark\\blst_benchmarks.exe --benchmark_filter=bench_g1_fft_trunc_middle_shape_detailed")
$lines.Add("out\\build\\x64-Release\\benchmark\\blst_benchmarks.exe --benchmark_filter=bench_gt_ifft_full_middle_shape_detailed")
$lines.Add("out\\build\\x64-Release\\benchmark\\blst_benchmarks.exe --benchmark_filter=bench_gt_ifft_trunc_middle_shape_detailed")
$lines.Add("out\\build\\x64-Release\\benchmark\\blst_benchmarks.exe --benchmark_filter=bench_middle_product_basic_no_trunc_detailed")
$lines.Add("out\\build\\x64-Release\\benchmark\\blst_benchmarks.exe --benchmark_filter=bench_middle_product_basic_ifft_trunc_detailed")
$lines.Add("out\\build\\x64-Release\\benchmark\\blst_benchmarks.exe --benchmark_filter=bench_middle_product_basic_pruned_detailed")
$lines.Add("out\\build\\x64-Release\\benchmark\\blst_benchmarks.exe --benchmark_filter=bench_middle_product_fft_detailed")
$lines.Add("out\\build\\x64-Release\\benchmark\\blst_benchmarks.exe --benchmark_filter=bench_middle_product_specialized_trunc_detailed")
$lines.Add("out\\build\\x64-Release\\benchmark\\blst_benchmarks.exe --benchmark_filter=bench_middle_product_one_ifft_detailed")
$lines.Add("\end{verbatim}")
$lines.Add("")
$lines.Add("\end{document}")

$outDir = Split-Path -Parent $OutTex
if (-not (Test-Path $outDir)) {
    New-Item -ItemType Directory -Path $outDir | Out-Null
}
$lines | Set-Content -Path $OutTex -Encoding ASCII

# Also emit compact CSVs for quick plotting.
$kernelCsv = @()
foreach ($r in ($g1FullRows | Sort-Object { [int]$_.n })) {
    $n = [int]$r.n
    $kernelCsv += [pscustomobject]@{
        n = $n
        m = [int]$r.m
        g1_fft_full_ms = [double]$r.real_time
        g1_fft_trunc_ms = [double]$g1TruncByN[$n].real_time
        gt_ifft_full_ms = [double]($gtFullRows | Where-Object { [int]$_.n -eq $n } | Select-Object -First 1).real_time
        gt_ifft_trunc_ms = [double]$gtTruncByN[$n].real_time
    }
}
$kernelCsv | Export-Csv -NoTypeInformation -Path "docs/planning/fft_kernel_sweep.csv"

$mpCsv = @()
foreach ($r in ($basicNoTruncRows | Sort-Object { [int]$_.n })) {
    $n = [int]$r.n
    $mpCsv += [pscustomobject]@{
        n = $n
        m = [int]$r.m
        basic_no_trunc_ms = Get-TotalMs $r
        basic_ifft_trunc_ms = Get-TotalMs $basicIfftTruncByN[$n]
        basic_pruned_ms = Get-TotalMs $basicPrunedByN[$n]
        specialized_ms = Get-TotalMs $specByN[$n]
        specialized_trunc_ms = Get-TotalMs $specTruncByN[$n]
        oneifft_cyclic_ms = Get-TotalMs $oneIfftByN[$n]
    }
}
$mpCsv | Export-Csv -NoTypeInformation -Path "docs/planning/middle_product_sweep.csv"

Write-Host "Wrote $OutTex"
Write-Host "Wrote docs/planning/fft_kernel_sweep.csv"
Write-Host "Wrote docs/planning/middle_product_sweep.csv"
