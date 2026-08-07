# File Review: `blst/scalar/fields/Fp.h`

| Field | Value |
| --- | --- |
| Disposition | in-scope production |
| Profiles | scalar, C, C++, Rust |
| Review state | in review |
| Exact review baseline | `535c4c6cdea0f6c31fafe01c305d7bb74015f610` |
| Remediation commits | `7c0e399`, `b6a51bb62fa4e52e8349e9d7a17d51617c212e24`, `61c8cd4`, `314c186`, `535c4c6` |
| Upstream lineage | `54e6e556:src/consts.c`, `src/consts.h`, `src/recip.c`, and addition-chain headers; no one-to-one file |
| Findings | `AUD-067`, `AUD-068`, `AUD-069`, `AUD-070`, `AUD-071`, `AUD-198`, `AUD-199`, `AUD-201`, `AUD-203` |

## Review result

Every line of the host Fp implementation, direct field callers, relevant C/Rust
boundaries, and focused tests was reviewed. Constants match the pinned upstream
literals. Fixed storage, unrolled addition chains, selected native calls, and
public-exponent control flow are accepted hot-path patterns. The host file has
no allocation or ordinary exception source.

Zero-bit exponentiation now returns one. Repeated square/multiply defines count
zero, protects exact multiplier aliasing, and uses full-width generated x86
counters. The C square-root wrapper is failure-atomic without changing the
candidate-producing core used by Fp2. Focused behavior, deterministic assembly
generation, and optimized native code-generation evidence pass.

The Montgomery identity is now one real `Fp12` object whose Fp and Fp2 views
are live nested subobjects. Reciprocal uses one 96-byte `UFp` limb array without
union lifetime switching. Its optimized native stack size and selected call
sequence are unchanged.

## Residual risk

`AUD-070` is validated. `AUD-069` now has a fixed zero-copy public/internal
type model, with
focused runtime, Rust, and Apple validation still pending. Retained Rust
execution for the changed square-root failure output, independent big-integer
constants, and full sanitizer/oracle evidence remain follow-ups.

## Generated-addition-chain follow-up

The `AUD-181` through `AUD-184` packet machine-checks the complete p-2 and
(p-3)/4 schedules, moves them into an Fp-only detail namespace, replaces their
initial byte copy with typed assignment, and pins their original BLST commits.
The native fixed-chain contract now states the zero-to-zero behavior required
by branchless infinity normalization; `recip_sqrt(0)` succeeds with output
zero. The public C fixed inverse retains its narrower nonzero precondition.

Scalar MSVC Debug and Release focused tests pass 3/3. Optimized reciprocal and
reciprocal-square-root frames and arithmetic call counts are unchanged; the
typed assignment removes one fixed pointer-identity branch from each entry.
The follow-up implementation is recorded at
`529a79377b6ef848ec90fce35e7580a9f451fc4f`. Generic scratch erasure was
explicitly declined and is not part of the constant-time contract.

## Wave 2 exact-contract follow-up

The complete file was reread at `535c4c6`. Public-exponent and public-count
control flow is now explicit, REDC states and Debug-checks its `a < p*R` bound,
and extreme public shifts use logarithmic work without changing the selected
kernel through count 256. An optimized count-three probe remains a direct
`lshift_mod_384` call. Retained Debug tests exercise both sides of REDC's upper
bound. The result-verification/fallback branch in `reciprocal` is intentional:
upstream added it for reported almost-inverse corner cases. Optimized
MSVC/Clang/GCC probes retain that branch, while `flt_reciprocal` calls only the
fixed chain. `AUD-201` is validated with the Euclidean operation documented as
variable-time rather than promoted to a stronger contract.
