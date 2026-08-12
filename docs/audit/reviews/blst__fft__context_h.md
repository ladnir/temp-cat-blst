# File Review: `blst/fft/context.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/fft/context.h` |
| Roadmap phase / queue | Phase 17 FFT support-header packet |
| Reviewed commit | `0fc26d2617f74b6ef82345fee4f4400979b1534d` |
| Reviewer | Codex audit team |
| Date | 2026-08-11 |
| Disposition | `in_scope_production` |
| Subsystem | Compile-time scalar-field, group, and pairing coefficient operations for generic FFTs |
| Profiles | scalar and selected AVX native C++; Fr C/Rust adapters |
| Generation/provenance | Locally composed handwritten C++; Apache-2.0 and Category Labs |
| Upstream references | `not_applicable: upstream BLST has no corresponding generic FFT context` |
| Related files | `blst/fft/ctFft.h`, `blst/fft/domain.h`, `blst/fft/middle_product.h`, selected `Fr`, G1, G2, and GT implementations, `test/Fft_api_tests.cpp`, `test/Fft_tests.cpp` |

## Purpose and trust boundary

`FrCoeffCtx` supplies Fr arithmetic to scalar transforms. `CoeffCtx` supplies
same-group G1/G2/GT addition, subtraction, and negation; public-Fr scalar action;
and G1-by-G2 pairing for polynomial products. Dispatch is compile-time through
the selected native types.

FFT shapes and roots are public. Coefficients and group values may be secret-
derived, but the group scalar action is deliberately variable-time and therefore
requires its Fr operand to be a public transform root or scale.

## Contracts and invariants

- Fr inputs are reduced Montgomery values. Inversion operands are nonzero.
- G1/G2 operands are on-curve prime-order subgroup elements; GT operands satisfy
  target-group membership. Release performs no membership validation here.
- `plus`, `minus`, and `neg` preserve the group domain.
- Group-by-Fr multiplication converts the public Montgomery Fr to canonical
  little-endian scalar bytes and calls `mult_vartime` at 256 bits.
- G1-by-G2 multiplication performs affine conversion, a Miller loop, and final
  exponentiation, returning GT.
- Results are returned by value; operations retain no references, allocate no
  storage, and contain no throwing construct. Arithmetic output temporaries are
  intentionally not value-initialized because each primitive overwrites them.

## Dependency and call-path map

The retained radix-two FFT calls `plus`, `minus`, `mult`, `inverse`, and
`from_u64` through a template parameter. Fr C/Rust transforms use `FrCoeffCtx`.
Native G1/G2/GT FFT tests and polynomial pairing products use the aliases
`G1CoeffCtx`, `G2CoeffCtx`, and `GTCoeffCtx`. Middle product owns more specialized
contexts for prepared AVX GT roots; this file does not add runtime selection or
global caching.

## Upstream/specification comparison

This context is local because upstream BLST has no generic FFT engine. Its group
contracts follow the audited CatBLST convention: ordinary G1/G2/GT operations
assume prime-order subgroup inputs unless explicitly documented otherwise.

## Checklist evidence

### 1. Identity and purpose

The complete header, history, aliases, selected-backend types, direct callers,
tests, and neighboring specialized contexts were reviewed.

### 2. Contracts and invariants

Representation, subgroup, public-scalar timing, nonzero inversion, result,
allocation, exception, and ownership contracts are documented in source.

### 3. Upstream and specification lineage

The local implementation now records Apache-2.0/Category Labs provenance and
the absence of an upstream BLST counterpart.

### 4. Correctness

Fr operations map directly to Fr primitives. Group subtraction negates the
right operand then uses complete group addition. Root action matches direct
selected-backend `mult_vartime`. Pairing the two generators returns the compiled
GT generator in scalar and portable-AVX-selected Debug profiles.

### 5. Memory and lifetime safety

All objects are fixed-size values. Canonical Fr bytes are read through the
permitted byte representation. No pointer escapes and no heap storage exists.

### 6. Undefined and implementation-defined behavior

The scalar byte view uses CatBLST's audited 64-bit little-endian Fr layout and
`byte` character-type alias. Unsupported template combinations fail the explicit
compile-time assertion instead of silently selecting a different operation.

### 7. Side channels and cryptographic hygiene

Fr arithmetic and complete group addition retain their established timing
contracts. Group scalar multiplication is explicitly variable-time and permitted
only because FFT roots/scales are public. No coefficient-dependent dispatch was
introduced.

### 8. Randomness and secret lifecycle

N/A: no randomness or persistent storage. The canonical scalar temporary holds
only a public root and therefore requires no secure erasure.

### 9. API and ABI

The types are native header templates/aliases and export no binary symbols.
The C/Rust public FFT surface reaches only `FrCoeffCtx` through its reviewed C
adapter; native group FFT contracts are documented here.

### 10. Concurrency

Every operation is stateless and reentrant. There is no cache, thread-local
selection, thread pool, or retained caller storage.

### 11. Resource exhaustion and denial of service

N/A: individual operations are fixed-size and allocation-free. Transform-size
resource policy belongs to `ctFft.h` and `domain.h`.

### 12. Portability and build behavior

The header names its standard-library dependencies and compiles when included
first by the focused test. Scalar-default and portable-AVX-selected MSVC Debug
profiles instantiate every retained operation family.

### 13. Performance-sensitive structure

Compile-time dispatch, direct primitive calls, canonical stack bytes, and
uninitialized fully overwritten outputs match the original code-generation
shape. The ambiguous equal-type branch was removed; no runtime check, virtual
call, allocation, type erasure, or extra value initialization was added.

### 14. Tests and evidence

Direct selected-backend checks cover G1, G2, and GT public-Fr multiplication and
generator pairing. Existing scalar group FFT round trips and convolution tests
cover addition/subtraction composition.

### Specialized additions

The generic name `CoeffCtx` is retained to avoid unnecessary source churn, but
its legal operand domains are now explicit. Specialized prepared AVX GT root
handling remains owned by `blst/avx/fft/gt_context.h`.

## Findings

`AUD-088` and `AUD-305`.

## Validation

- MSVC 19.50 scalar-default Debug focused FFT/context tests: passed.
- MSVC 19.50 portable-AVX-selected Debug focused domain/context tests: 4/4
  passed.
- Final exact source compiles in both selected profiles after removal of
  unnecessary temporary initialization.
- No benchmark, CUDA, full suite, sanitizer, native AVX-512, remote, or Apple
  run was performed.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Earlier Phase 12 database | `AUD-088` timing-contract context | Cannot prove subgroup mathematics or constant-time behavior |
| Clang Static Analyzer / clang-tidy | Earlier retained profiles | No mapped diagnostic | Exact commit validated with MSVC |
| Clippy | N/A | Rust adapter reviewed separately | Native C++ template header |
| ASan/UBSan | Earlier Phase 17 focused run | No mapped context failure | Not rerun for contract/type-only refactor |
| Fuzzing | N/A | N/A | Fixed-size adapter; arithmetic owners have their own differential evidence |

## Residual risk and decision

The context relies on caller-established subgroup membership and public roots;
Release validation would violate the intended narrow, high-performance contract.
Native AVX-512 hardware and remaining milestone evidence are pending. The file
is complete for line-by-line hand review.
