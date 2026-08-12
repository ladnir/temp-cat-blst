# File Review: `bindings/blst_avx.h`

## Record

| Field | Value |
| --- | --- |
| Roadmap phase | Phase 18 expert AVX public surface |
| Reviewed commit | `e6163cecba9a9f66a8e6ec7aed58cbf2f0ba1672` |
| Date | 2026-08-12 |
| Disposition | in-scope production, complete |
| Profiles | C99, C++, native/portable AVX-512; CUDA and frozen BTx/BTE excluded |
| Related findings | `AUD-091`–`AUD-095`, `AUD-160`, `AUD-161`, `AUD-276`, `AUD-316`–`AUD-322` |

## Purpose and trust boundary

This header is the expert C ABI for opaque AVX point and Fp12 carriers,
representation conversion, basic point operations, subgroup multiplication and
MSM, encoding, membership, final exponentiation, and pairing. Callers own every
object and scratch buffer. Opaque bytes are an ABI carrier, not a wire format or
an invitation to construct field values manually.

## Contracts and review results

- All 77 declarations were inventoried against definitions and handwritten Rust
  declarations. Rust now has exact 77/77 parity.
- Basic point operations document their on-curve domain. Optimized scalar and
  multipoint operations require prime-order subgroup elements and public
  scalars where variable-time behavior applies; Release does not validate those
  preconditions.
- Generic Fp12 multiply and square accept arbitrary well-formed Fp12 carriers,
  do not require or establish GT, and permit exact output aliasing. GT MSM is a
  separate narrow operation.
- Final exponentiation accepts a valid nonzero Miller-loop Fp12 value and
  produces GT. Its C representation does not make MillerResult and GT the same
  C++ type.
- Same-type arithmetic permits exact aliasing only where stated. Partial overlap
  and differently typed conversion overlap are forbidden.
- Scratch is caller-owned, sized by the matching helpers, and exclusively live
  for the call. Zero can mean empty work or invalid/unavailable sizing as
  documented by the operation.
- The six G1 GLV scratch/window hooks are retained. BTX-DKG production uses the
  automatic family and its tuning benchmarks use the explicit-window family.
- No entry allocates or permits an exception to cross C. Callback exceptions are
  translated to failure by the implementation.
- The ABI sizes and 64-byte alignment remain unchanged. No runtime dispatch,
  dynamic polymorphism, or hidden conversion layer was introduced.

## Evidence and residual risk

Strict MSVC C11 and C++20 `/W4 /WX /Zs` parses pass. Portable-AVX Debug passes
the generic-Fp12/final-exponent and retained-GLV-hook regressions. Rust AVX tests
pass 7/7. MatrixNoK dispatch and 21 accepted correspondence rows/61 hashes
remain valid. Native AVX-512 runtime, Apple, sanitizer/fuzz, installed-SDK, and
broader cross-compiler evidence remain phase-wide gates. No benchmark was run
because the only arithmetic optimization selects the existing dedicated field
square kernel; native performance confirmation remains appropriate at the
milestone gate.
