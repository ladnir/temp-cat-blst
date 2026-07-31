# Phase 7B C2x Public-Header Compatibility Evidence

## Record

- Date: 2026-07-31
- Implementation commit: `d9aed1e`
- State: validated and committed
- Scope: upstream `f262a6e`
- Benchmarks: excluded and not run
- CUDA: experimental and not exercised
- BTx/BTE: frozen and untouched

## Decision

Translate the upstream fix into `bindings/blst_types.h`, where the fork's split
public headers centralize C-compatible types.

For ordinary C99 through pre-final C23 modes, the header now includes
`<stdbool.h>` and lets the compiler's standard library select the correct
definition. The Rust-bindgen preprocessing mode retains an explicit
`bool -> _Bool` macro because that generator depends on it. The legacy
pre-C99 fallback remains `bool -> int`, and C++ and final C23 modes remain
untouched.

This removes the fork's unconditional `_Bool` macro substitution in draft-C2x
mode, which was the ambiguity fixed upstream, without restoring the removed Go
or Zig branches.

## Validation

Following `VALIDATION_POLICY.md`, validation was limited to the public-header
surface:

| Surface | Configuration | Result |
| --- | --- | --- |
| C header | GCC 13.3, `-std=c89` legacy fallback | Passed |
| C header | GCC 13.3, `-std=c99`, `c11`, and `c17` with warnings as errors | Passed |
| Draft C2x header | GCC 13.3, `-std=c2x` with warnings as errors | Passed |
| Rust-bindgen preprocessing contract | GCC 13.3, C2x with `__BLST_RUST_BINDGEN__` | Passed |
| Native C++ header | G++ 13.3, C++23 with warnings as errors | Passed |
| Formatting/integrity | Git whitespace and audit manifest checks | Passed |

No product library, AVX-512, sanitizer, Rust test suite, CUDA, or benchmark run
was warranted because the change is entirely compile-time public-header
selection and does not alter declarations, layout, symbols, or runtime code.

## Residual validation

- The available GCC identifies the ratification-era mode as `c2x`; a compiler
  exposing final `__STDC_VERSION__ == 202311L` was not available locally.
- Clang and the bindgen executable were not available on this host. Their next
  planned CI or Rust-generation gate should repeat the same header check.
- MSVC has no final C23 mode in the available toolchain; its existing legacy C
  fallback is unchanged.
