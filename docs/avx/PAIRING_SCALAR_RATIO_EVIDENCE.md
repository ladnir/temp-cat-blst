# MatrixNoK Pairing-to-Scalar Evidence

## Result

On Peach's AMD Ryzen 9 7950X, the selected production MatrixNoK path is 3.59x
faster than scalar CatBLST for the Miller loop and 3.53x faster for a complete
pairing. The selected complete-pairing wrapper is within 0.55% of calling the
MatrixNoK engine directly, so the result does not expose material integration
overhead.

The result is below VROOM's approximately 4x Sapphire Rapids result because
the scalar timings are similar while Peach executes the vector path more
slowly. It does not indicate that CatBLST failed to preserve VROOM's selected
engine performance: the selected CatBLST pairing remains approximately equal
to the pinned VROOM pairing previously measured on the same Peach host.

## Configuration

- Production source commit: `d314a3a`.
- Host: AMD Ryzen 9 7950X, native AVX-512IFMA.
- Compiler: Clang 21.1.8.
- Release flags: `-O3 -march=native -DNDEBUG`.
- Scalar configuration: assembly and ADX enabled.
- Selected vector configuration: native AVX-512, MatrixNoK default, portable
  AVX disabled.
- CUDA, sanitizers, tests, parallel MSM, and portable AVX were disabled.
- One process ran at a time, pinned to CPU 2.
- Nine repetitions used a 0.2-second minimum case time and report medians.
- CPU scaling remained enabled, so these are engineering-selection results,
  not final release claims.

## Selected public-path comparison

| Operation | Scalar | Selected MatrixNoK | Speedup |
| --- | ---: | ---: | ---: |
| Miller loop | 208.281 us | 58.008 us | 3.590x |
| Complete pairing | 492.038 us | 139.388 us | 3.530x |

For context, VROOM's camera-ready Sapphire Rapids artifact reports 215.023 us
versus 49.888 us for Miller (4.310x) and 500.035 us versus 128.146 us for a
complete pairing (3.902x). Peach's scalar result is therefore not materially
faster than the paper baseline. The difference is concentrated in the native
vector result, consistent with Peach's split 512-bit execution.

## Integration-overhead diagnostic

A second serial process compared the selected public path with the same
MatrixNoK templates called directly in the same executable:

| Operation | Selected public path | Direct engine | Difference |
| --- | ---: | ---: | ---: |
| Miller loop | 58.086 us | 55.893 us | +3.92% |
| Complete pairing | 138.214 us | 137.454 us | +0.55% |

The direct Miller benchmark constructs its `Miller` owner before entering the
timed loop, while the public operation constructs it per invocation. That
explains why the Miller-only delta is not an equivalent wrapper-overhead
measurement. Complete pairing includes corresponding setup on both paths and
is the relevant integration gate.

The previously recorded pinned-VROOM complete pairing on Peach was 139.728 us.
The current selected CatBLST measurements of 138.2--139.4 us remain at parity
or slightly faster. There is no measured VROOM implementation gap to recover
on Peach; reaching VROOM's approximately 4x scalar ratio requires the stronger
512-bit execution of its Sapphire Rapids reference host.

## Benchmark infrastructure

`avx_pairing_ratio_benchmarks` is deliberately independent of the broad
benchmark executable and legacy test helpers. It measures only production
scalar, selected AVX, and direct selected-engine pairing paths. This keeps the
gate usable while the larger benchmark tree is updated after the core type
redesign.
