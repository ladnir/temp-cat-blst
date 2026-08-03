# File Review: `blst/scalar/fields/Fp.h`

| Field | Value |
| --- | --- |
| Disposition | in-scope production |
| Profiles | scalar, C, C++, Rust |
| Review state | in review |
| Exact review baseline | `c3096482cdab88edc4f05e0c6583c7896fd75bab` |
| Remediation commits | `7c0e399`, `b6a51bb62fa4e52e8349e9d7a17d51617c212e24`, `61c8cd4` |
| Upstream lineage | `54e6e556:src/consts.c`, `src/consts.h`, `src/recip.c`, and addition-chain headers; no one-to-one file |
| Findings | `AUD-067`, `AUD-068`, `AUD-069`, `AUD-070`, `AUD-071` |

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
