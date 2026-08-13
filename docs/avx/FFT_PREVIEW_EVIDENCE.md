# MatrixNoK FFT Preview Evidence

## Decision

The pinned MatrixNoK field tower is worth porting through G1 and GT. In a
same-process native AVX-512IFMA preview, its forward G1 FFT has 9.3% to 11.0%
lower latency than CatBLST's incumbent 50-bit AVX path, equivalent to a
1.10x to 1.12x speedup. Its inverse GT FFT has 15.8% to 16.9% lower latency,
equivalent to a 1.19x to 1.20x speedup.

These historical preview results did not change dispatch at the time. MatrixNoK
now owns the production AVX field, G1, and GT paths.

## Construction

- CatBLST production arithmetic state:
  `e260936f85f07ac82fa363433a754519f0e7a6d1`.
- Benchmark-only preview adapter:
  `71bfbfdefffbdd630bbf83bc774a8f73938cb8c1`.
- Corrected G1 adapter and non-circular validity gate:
  `edaf5462fd95e40bb6653f138b0562f9f85130a8`.
- VROOM source: exact commit
  `26bd9aee2168f167583cf096add355547a7b56d4`.
- Host: AMD Ryzen 9 7950X, native AVX-512IFMA, one process pinned to CPU 2.
- Compiler: Ubuntu Clang 21.1.8, Release, `-O3 -march=native`.
- CUDA, BTx/BTE, portable-AVX Release, sanitizers, and broad tests were not
  built or run.

CatBLST and VROOM are compiled in separate translation units of the existing
`avx_group_ops_benchmarks` executable. This is required because VROOM's
headers define process-global parameter objects. A benchmark-only bridge
supplies canonical Fr roots, compiled CatBLST GLV4 schedules, canonical G1
generator coordinates, scalar GT generator coefficients, and VROOM's
inherited `blst_div_by_zz` spelling.

Both paths execute the same radix-2 butterfly schedule. Root preparation,
input construction, forward preparation for the inverse case, and per-run
input reset are outside the timed region.

## Correctness gates

The GT adapter imports CatBLST's GT generator as 12 canonical base-field
coefficients. Before timing, it performs a complete 256-element forward and
inverse FFT round trip with the pinned VROOM Fp12 implementation and the same
public-root GLV4 schedules used by CatBLST. The gate passes.

The G1 adapter is deliberately fail-closed. The original preview copied the
G1 generator as a decimal literal and mistyped one digit in `x`, producing an
off-curve input. Its VROOM-versus-VROOM multiplier check did not detect the
error because both sides consumed the same invalid point. The resulting
addition failures were adapter failures, not VROOM curve failures.

The corrected adapter imports CatBLST's canonical generator coordinates and
checks the homogeneous curve equation independently with `BigInt` before any
VROOM curve operation is accepted. It uses VROOM's EC-qualified MatrixNoK
ring and production `add_point` kernel. Before timing, it passes:

- generator and doubled-generator curve-equation checks;
- the two-element identities `P + [2]P = [3]P` and
  `P - [2]P = -P`;
- public GLV agreement for zero, one, two, representative roots, and the
  inverse-of-two scale; and
- a complete 256-element forward/inverse G1 FFT round trip.

This incident demonstrates that same-implementation differential checks are
not sufficient input validation. Candidate curve tests must retain an
independent scalar or mathematical oracle.

## Native G1 results

The first three sizes use three same-process repetitions and report the mean.
The 16,384 case uses one repetition. Google Benchmark warned that CPU scaling
was enabled, so these are port-selection results rather than final release
claims.

| Elements | CatBLST incumbent | MatrixNoK preview | Latency reduction | Speedup |
| ---: | ---: | ---: | ---: | ---: |
| 256 | 20.5 ms | 18.6 ms | 9.3% | 1.102x |
| 1,024 | 111 ms | 100 ms | 9.9% | 1.110x |
| 4,096 | 564 ms | 504 ms | 10.6% | 1.119x |
| 16,384 | 2,735 ms | 2,433 ms | 11.0% | 1.124x |

The repeated G1 cases reported coefficient of variation at or below 0.02%.

## Native GT results

The first three sizes use three randomized-interleaved repetitions and report
the median. The 16,384 case uses one repetition to keep the approved remote
job bounded. Google Benchmark warned that CPU scaling was enabled, so the
result is suitable for a port/no-port decision but not a final release claim.

| Elements | CatBLST incumbent | MatrixNoK preview | Latency reduction | Speedup |
| ---: | ---: | ---: | ---: | ---: |
| 256 | 36.1 ms | 30.0 ms | 16.9% | 1.203x |
| 1,024 | 190 ms | 158 ms | 16.8% | 1.203x |
| 4,096 | 941 ms | 784 ms | 16.7% | 1.200x |
| 16,384 | 4,469 ms | 3,761 ms | 15.8% | 1.188x |

The three repeated MatrixNoK cases had reported coefficient of variation at
or below 0.02%; the incumbent cases were also at or below 0.02%.

## Historical interpretation and completed gate

GT FFT is dominated by Fp12 multiplication, cyclotomic squaring, Frobenius,
and public-root exponentiation. The consistent end-to-end gain is stronger
evidence than the earlier isolated Fp/Fp2 timings and justifies completing the
Fp12/GT candidate port.

Later packets ported the audited CatBLST curve and Fp12/GT operations onto
MatrixNoK, added independent scalar differentials, and selected the resulting
production path. VROOM established the performance target; CatBLST retained
its reviewed point contracts rather than importing the research API wholesale.
