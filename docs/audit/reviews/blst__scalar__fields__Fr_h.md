# File Review: `blst/scalar/fields/Fr.h`

| Field | Value |
| --- | --- |
| Disposition | in-scope production |
| Profiles | scalar, C, C++, Rust |
| Review state | in review |
| Exact review baseline | `314c186149b9cc95d85a0c54e2dd651f137e4236` |
| Remediation commits | `b6a51bb62fa4e52e8349e9d7a17d51617c212e24`, `61c8cd4`, `314c186` |
| Upstream lineage | `54e6e556:src/consts.c`, `src/exports.c`, `src/recip.c`, and `src/vect.c`; no one-to-one file |
| Findings | `AUD-067`, `AUD-069`, `AUD-070`, `AUD-072`, `AUD-198`, `AUD-199`, `AUD-201`, `AUD-203` |

## Review result

Every line of the host Fr implementation, its direct callers, supported C/Rust
boundaries, and relevant tests was reviewed. Constants and selected arithmetic
structure match the pinned upstream lineage. The host implementation allocates
nothing and contains no ordinary exception source. Public-exponent branching,
fixed arrays, forced-inline arithmetic, and selected assembly calls are
accepted performance-sensitive patterns.

Zero-bit exponentiation now returns Montgomery one before reading the exponent.
Secret-boundary temporaries are explicitly cleared, dormant 32-bit conversion
code is removed, and representation documentation/provenance is corrected.
Focused behavior, strict compilation, and optimized wipe/identity code
generation pass.

Canonical scalar conversion now uses the documented little-endian platform
policy, exact typed in-place Montgomery conversion, and fixed-size copies into
live Fr temporaries. Inactive union reads, whole-object casts, and redundant
alignment branches are removed.

## Residual risk

`AUD-070` is validated. `AUD-069` now has a fixed zero-copy
public/internal type model, with focused runtime, Rust, and Apple validation
still pending. Butterfly nonalias contracts, inverse of
zero, and full sanitizer/oracle evidence remain Phase 14 follow-ups.

## Wave 2 exact-contract follow-up

The complete file was reread at `314c186`. Exponent pointer/top-bit and public
timing requirements are explicit. REDC now states `a < r*R` and checks the
upper-half bound only in Debug. Counts through 256 retain the original shift
backend; larger public counts use logarithmic modular exponentiation, including
`SIZE_MAX`. Constant-count optimized probes remain direct calls to the original
backend. Cross-platform and independent-oracle evidence remains phase work.
