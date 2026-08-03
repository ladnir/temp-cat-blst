# File Review: `blst/scalar/fields/Fp.h`

| Field | Value |
| --- | --- |
| Disposition | in-scope production |
| Profiles | scalar, C, C++, Rust |
| Review state | in review |
| Exact review baseline | `c3096482cdab88edc4f05e0c6583c7896fd75bab` |
| Remediation commits | `7c0e399`, `b6a51bb62fa4e52e8349e9d7a17d51617c212e24` |
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

## Residual risk

The file remains in review for `AUD-070` typed constant/reciprocal lifetime
design and the wider `AUD-069` C storage boundary. Retained Rust execution for
the changed square-root failure output, Apple execution, independent big-
integer constants, and full sanitizer/oracle evidence remain follow-ups.
