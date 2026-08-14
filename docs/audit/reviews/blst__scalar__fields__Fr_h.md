# File Review: `blst/scalar/fields/Fr.h`

| Field | Value |
| --- | --- |
| Disposition | in-scope production |
| Profiles | scalar, C, C++, Rust |
| Review state | in review |
| Exact review baseline | `535c4c6cdea0f6c31fafe01c305d7bb74015f610` |
| Remediation commits | `b6a51bb62fa4e52e8349e9d7a17d51617c212e24`, `61c8cd4`, `314c186`, `535c4c6`, `1e286e2` |
| Upstream lineage | `54e6e556:src/consts.c`, `src/exports.c`, `src/recip.c`, and `src/vect.c`; no one-to-one file |
| Findings | `AUD-067`, `AUD-069`, `AUD-070`, `AUD-072`, `AUD-198`, `AUD-199`, `AUD-201`, `AUD-203`, `AUD-222` |

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

`AUD-069` and `AUD-070` are validated. Exact C/C++ type assertions, focused
boundary behavior, and supported native C++/Rust run `31665323681` cover the
zero-copy public/internal representation on Apple Silicon, Linux, and Windows.
The later assembly-layout and independent-oracle additions at `a29895f` and
`83c11ba` still require Apple execution. Broad sanitizer evidence remains a
phase gate.

## Wave 2 exact-contract follow-up

The complete file was reread at `535c4c6`. Exponent pointer/top-bit and public
timing requirements are explicit. REDC now states `a < r*R` and checks the
upper-half bound only in Debug. Counts through 256 retain the original shift
backend; larger public counts use logarithmic modular exponentiation, including
`SIZE_MAX`. Constant-count optimized probes remain direct calls to the original
backend. Retained Debug tests cover the exact REDC upper bound. The selected Fr
inverse reaches fixed-count generated assembly on x86 and fixed-iteration
portable source; optimized MSVC/Clang/GCC wrappers add no value-dependent
control flow. Independent-oracle and broader runtime evidence remains phase
work.

`AUD-222` distinguishes two reduction contracts that the original call graph
conflated. Public `Fr::redc` remains restricted to `T < rR` and publishes a
reduced field element. The private almost-inverse is bounded by the
left-aligned modulus `2r`, may produce a first-stage value below `2r`, and is
canonicalized by the immediately following `to_mont`. Debug now checks that
private bound without adding Release instructions.

## Assembly-layout follow-up

At `a29895f`, Fr is required at compile time to remain a 32-byte standard-layout,
trivially copyable object whose limb array starts at offset zero. The independent
Montgomery oracle covers reduced operands, the sparse full-width second operand,
REDC, from-Montgomery, and exact aliases. `AUD-371` remains fixed pending the
supported Apple compile.

## Batch D inverse follow-up

At `83c11ba`, the private bound comment records both `temp * inp = R^2` and
`temp < 2rR`. The existing direct REDC plus `to_mont` pipeline is unchanged.
Independent selected and portable x86 oracles validate zero, the bound, and
the congruence. Apple execution remains pending under `AUD-373`.
