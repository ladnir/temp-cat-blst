# File Review: `bindings/blst_avx.h`

## Record

| Field | Value |
| --- | --- |
| Roadmap phase | Phase 18 expert AVX public surface |
| Reviewed commit | `daa4eb3005fa4c2fef68fd23a84c34dd2518df7d` |
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

## Phase 16 MSM boundary addendum

Commit `5b08edc8f890ef1cbf4c6443162ec4da985879ee` documents the exact live-range
contract for affine G1/G2/GT MSM, expert G1 GLV, and projective multipoint MSM.
Required writable scratch and all input/output ranges are pairwise disjoint.
One-point projective exact output/input aliasing remains supported. Invalid
overlap returns `false` without changing output; entries allocate nothing and
do not throw across C. Focused portable-Debug tests validate exact and partial
overlap, failure atomicity, and one-point G1/G2 aliasing. `AUD-390` is
validated; native, sanitizer, and fuzz evidence remain phase-wide gates.

## Pairing result addendum

At `daa4eb30`, `blst_avx_pairing` remains a subgroup-only expert operation.
Its boolean reports null arguments or unavailable AVX support; it does not
report curve or subgroup validation. No implementation or code generation
changed. This closes the AVX declaration portion of `AUD-097`.

## 2026-08-27 secret G1 MSM addendum

The public C surface now separates constant-time secret MSM from the existing
variable-time Pippenger families. `blst_avx_p1s_mult_ct` documents resident
projective subgroup points, reduced native `blst_fr` scalars, public count,
fixed 32/16/8 chunking, scratch size/alignment/disjointness, empty input, and
failure atomicity. `blst_avx_p1s_mult_ct_scratch_sizeof` reports the maximum
reusable scratch requirement. The API allocates nothing and makes no
non-subgroup or noncanonical-scalar promise.

## 2026-08-27 public native-scalar MSM addendum

The header adds two declarations, bringing the reviewed C surface to 79
functions. `blst_avx_p1s_mult_fr_vartime` explicitly restricts its resident
`blst_fr` inputs to public reduced scalars and documents subgroup ownership,
allocation freedom, scratch alignment and disjointness, empty input, and
failure atomicity. Its matching capacity query documents that the algorithm
may select shared-doubling GLV wNAF or Pippenger without exposing that policy as
a caller-visible representation change.

## 2026-09-04 alpha-pairing range addendum

The header adds a capacity-bounded prepared-G2 workspace constructor and one
contiguous pairing-range operation. The range documentation fixes the native
point provenance, nonoverlap, immutable-batch sharing, distinct-workspace,
allocation, public-timing, and failure-output contracts. These are kernel
operations only; no scheduler or phase state crosses the C boundary.

The reviewed C surface is now 100 functions. `DOC-028` maps both declarations
to their C++ implementation and safe Rust owner. Focused native-AVX Release and
ASan/UBSan tests pass on Peach.
