# BTX FFT and Middle-Product Baseline

## Purpose

This record is the Batch 0 baseline for the BTX FFT layout and range-scheduling
work. It fixes the comparison point before changing FFT schedule layout,
workspace ownership, or execution granularity. Later batches compare against
the medians below; this record is evidence, not an active implementation plan.

The measured implementation commit is
`4bd972107603daac6db6d8f6dadcb3f24fdee911`. The benchmark-only additions that
produce this report add the 1,024-point anchor, expose memory counters, and
build the FFT benchmark independently of unrelated benchmark translation
units. They do not change a CatBLST kernel.

## Environment and method

- Host: Peach, AMD Ryzen 9 7950X, 16 enabled cores, SMT disabled.
- OS: Linux 6.17.0-41-generic, x86-64.
- Compiler: Ubuntu Clang 21.1.8.
- Build: Release, `-O3 -march=native -DNDEBUG`, assembly enabled, native
  AVX-512 enabled, and `BLST_DEFAULT_AVX=ON`.
- Affinity: every benchmark process was pinned to CPU 4 with `taskset`.
- Isolation: no two benchmark processes ran concurrently.
- C++ sampling: five repetitions per case, with at least 100 ms selected by
  Google Benchmark. The tables report paired-run medians.
- Rust sampling: Criterion used a 500 ms warm-up and a two-second requested
  measurement interval. The table reports the center estimate and its
  confidence interval.

The host reported CPU frequency scaling as enabled. The five C++ repetitions
were nevertheless tightly grouped: the total-time coefficient of variation
was at most 0.08% for the scalar cases shown and at most 0.15% for the AVX
cases.

The C++ measurement command was:

```text
taskset -c 4 ./out/baseline-clang21/benchmark/fft_middle_product_benchmarks \
  --benchmark_filter='^bench_middle_product_(scalar|avx)_detailed/(32|64|128|256|512|1024)/real_time$' \
  --benchmark_min_time=0.1s --benchmark_repetitions=5 \
  --benchmark_report_aggregates_only=true --benchmark_out_format=json
```

The Rust boundary command used the same core and compiler configuration:

```text
CC=/usr/bin/clang-21 CXX=/usr/bin/clang++-21 \
CFLAGS='-O3 -march=native' CXXFLAGS='-O3 -march=native' \
RUSTFLAGS='-C target-cpu=native' \
taskset -c 4 cargo bench --manifest-path bindings/rust/Cargo.toml \
  --bench btx_native_middle_product --features avx
```

## Direct C++ middle-product baseline

All times are milliseconds. `n` is the logical middle-product size and the FFT
size is `2n` for these anchors. `Input` is identity fill plus input copy.
`Pointwise` includes the batch G1 affine conversion and all prepared pairings.
`Inverse GT` includes inverse-transform normalization. The current one-shot
implementation does not expose affine conversion or normalization as separate
timing regions; Batch 3 will expose those boundaries as independent phases.

### Native scalar

| n | Total | Input | Forward G1 | Pointwise | Inverse GT | Extract | Workspace bytes | FFT-plan bytes |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 32 | 56.644 | 0.000251 | 6.721 | 31.417 | 18.477 | 0.018 | 46,080 | 4,680 |
| 64 | 126.741 | 0.000462 | 17.220 | 62.798 | 46.652 | 0.036 | 92,160 | 10,408 |
| 128 | 280.895 | 0.000921 | 42.184 | 125.633 | 112.950 | 0.072 | 184,320 | 23,304 |
| 256 | 616.752 | 0.001763 | 99.657 | 251.300 | 265.564 | 0.145 | 368,640 | 52,072 |
| 512 | 1,343.376 | 0.003656 | 230.070 | 502.524 | 610.358 | 0.331 | 737,280 | 115,656 |
| 1,024 | 2,908.300 | 0.007464 | 522.128 | 1,005.134 | 1,380.000 | 0.653 | 1,474,560 | 255,016 |

### Native AVX-512

| n | Total | Input | Forward G1 | Pointwise | Inverse GT | Extract | Workspace bytes | FFT-plan bytes |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 32 | 13.929 | 0.000453 | 2.545 | 7.056 | 4.326 | 0.001 | 147,456 | 4,680 |
| 64 | 31.498 | 0.000942 | 6.511 | 14.108 | 10.873 | 0.002 | 294,912 | 10,408 |
| 128 | 70.398 | 0.002029 | 15.914 | 28.216 | 26.250 | 0.004 | 589,824 | 23,304 |
| 256 | 156.354 | 0.007283 | 37.646 | 57.020 | 61.607 | 0.017 | 1,179,648 | 52,072 |
| 512 | 343.356 | 0.019226 | 86.963 | 114.645 | 141.602 | 0.034 | 2,359,296 | 115,656 |
| 1,024 | 746.993 | 0.041328 | 197.217 | 229.836 | 319.551 | 0.060 | 4,718,592 | 255,016 |

The workspace figures are the retained heap payloads of the G1, affine-G1,
affine-prefix, and GT vectors after one evaluation. They exclude vector-object
headers and allocator metadata. The FFT-plan figures include vector payloads
for stage descriptors and nested index arrays; they exclude vector-object
headers and immutable roots and G2 precomputation. These definitions remain
fixed for the Batch 1 layout comparison.

## Rust native-AVX boundary check

The Rust benchmark calls the same native AVX precompute and evaluation path.
Its results agree with the direct C++ medians closely enough that no binding
overhead is visible at this scale.

| n | Criterion estimate | 95% interval |
| ---: | ---: | ---: |
| 32 | 13.957 ms | 13.955--13.959 ms |
| 64 | 31.554 ms | 31.552--31.555 ms |
| 128 | 70.532 ms | 70.528--70.540 ms |
| 256 | 156.60 ms | 156.57--156.63 ms |
| 512 | 344.12 ms | 344.08--344.17 ms |
| 1,024 | 749.04 ms | 748.94--749.17 ms |

## Validation

The Peach Clang 21 Release build completed the full CatBLST test executable:
336 tests passed and four tests were skipped because the optional shared BTX
vector path was not configured. No test failed.

Batch 1 must preserve scalar and AVX canonical outputs and compare sequential
middle-product medians against this record. A repeatable regression greater
than 2% fails the Batch 1 gate.
