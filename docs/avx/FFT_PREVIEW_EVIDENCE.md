# MatrixNoK FFT Preview Evidence

## Decision

The pinned MatrixNoK field tower is worth porting through GT. In a
same-process native AVX-512IFMA preview, its inverse GT FFT has 15.8% to 16.9%
lower latency than CatBLST's incumbent 50-bit AVX path, equivalent to a
1.19x to 1.20x speedup over the measured sizes.

No G1 performance claim is accepted. The pinned VROOM curve layer fails a
two-point G1 FFT base case under MatrixNoK, so CatBLST must first place its
audited curve formulas and bound handling over the candidate ring.

Neither result changes production dispatch.

## Construction

- CatBLST production arithmetic state:
  `e260936f85f07ac82fa363433a754519f0e7a6d1`.
- Benchmark-only preview adapter:
  `71bfbfdefffbdd630bbf83bc774a8f73938cb8c1`.
- VROOM source: exact commit
  `26bd9aee2168f167583cf096add355547a7b56d4`.
- Host: AMD Ryzen 9 7950X, native AVX-512IFMA, one process pinned to CPU 2.
- Compiler: Ubuntu Clang 21.1.8, Release, `-O3 -march=native`.
- CUDA, BTx/BTE, portable-AVX Release, sanitizers, and broad tests were not
  built or run.

CatBLST and VROOM are compiled in separate translation units of the existing
`avx_group_ops_benchmarks` executable. This is required because VROOM's
headers define process-global parameter objects. A benchmark-only bridge
supplies canonical Fr roots, compiled CatBLST GLV4 schedules, the scalar GT
generator coefficients, and VROOM's inherited `blst_div_by_zz` spelling.

Both paths execute the same radix-2 butterfly schedule. Root preparation,
input construction, forward preparation for the inverse case, and per-run
input reset are outside the timed region.

## Correctness gates

The GT adapter imports CatBLST's GT generator as 12 canonical base-field
coefficients. Before timing, it performs a complete 256-element forward and
inverse FFT round trip with the pinned VROOM Fp12 implementation and the same
public-root GLV4 schedules used by CatBLST. The gate passes.

The G1 adapter is deliberately fail-closed. Its public GLV executor agrees
with VROOM's independent BigInt overload for zero, one, two, representative
roots, and the inverse-of-two scale. The curve layer does not pass the
two-element FFT base case:

- An initial setup pass incorrectly used narrow `add_point` for `P + P`; the
  adapter was corrected to use VROOM's dedicated `double_point` operation.
- With that contract error removed, both `add_point` and `add_point_small`
  still fail the valid distinct-point subtraction identity
  `P - [2]P = -P`.

Consequently the G1 benchmark reports a setup error and produces no timing.
This is an upstream research implementation issue or an unresolved
MatrixNoK-bound integration issue, not evidence against the base-field
candidate. It is also why CatBLST's audited curve layer must be retained
rather than importing VROOM's curve file wholesale.

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

## Interpretation and next gate

GT FFT is dominated by Fp12 multiplication, cyclotomic squaring, Frobenius,
and public-root exponentiation. The consistent end-to-end gain is stronger
evidence than the earlier isolated Fp/Fp2 timings and justifies completing the
Fp12/GT candidate port.

The next packet should port the audited CatBLST Fp12 and GT operations onto
MatrixNoK, add scalar differential tests for every retained operation, and
repeat this benchmark against production-shaped code. G1 should follow as a
separate packet: first port CatBLST's reviewed point formulas and explicit
`add` versus `add_or_double` contracts, then pass focused curve differentials
and the two-point/256-point FFT gates before any performance run.
