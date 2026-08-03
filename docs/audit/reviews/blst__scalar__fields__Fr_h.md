# File Review: `blst/scalar/fields/Fr.h`

| Field | Value |
| --- | --- |
| Disposition | in-scope production |
| Profiles | scalar, C, C++, Rust |
| Review state | in review |
| Exact review baseline | `c3096482cdab88edc4f05e0c6583c7896fd75bab` |
| Remediation commit | `b6a51bb62fa4e52e8349e9d7a17d51617c212e24` |
| Upstream lineage | `54e6e556:src/consts.c`, `src/exports.c`, `src/recip.c`, and `src/vect.c`; no one-to-one file |
| Findings | `AUD-067`, `AUD-069`, `AUD-070`, `AUD-072` |

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

## Residual risk

The file remains in review because `AUD-070` covers inherited inactive-union
and in-place object-lifetime patterns, and `AUD-069` covers the wider public C
storage boundary. Butterfly nonalias contracts, big-endian exclusion, inverse
of zero, and full sanitizer/oracle evidence remain Phase 14 follow-ups.
