# File Review: `blst/avx/curve/multi_scalar.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/avx/curve/multi_scalar.h` |
| Roadmap phase / queue | Phase 16 AVX curve and MSM closure |
| Reviewed commit | `4c58255b4ea345fc11d367c98ebd4576a4773de2` |
| Reviewer | Codex audit team |
| Date | 2026-08-13 |
| Disposition | `in_scope_production` |
| Subsystem | AVX G1/G2/GT variable-time MSM |
| Profiles | C boundary, native AVX-512, portable developer Debug |
| Generation/provenance | handwritten CatBLST AVX MSM over selected curve types |
| Upstream references | not applicable: CatBLST AVX MSM |
| Related files | `blst/avx_c_api.cpp`, `bindings/blst_avx.h`, `Avx_msm_tests.cpp`, `Avx_c_api_contract_tests.cpp` |

## Purpose and trust boundary

The file owns variable-time public-scalar Pippenger kernels, bucket schedules,
GLV/GLS/GLV4 preparation, and scratch-size calculations. Ordinary G1/G2 inputs
must already lie in their prime-order subgroups. GT inputs must be valid GT
elements. The algorithms perform no Release membership validation.

Retained executors consume explicit caller-owned bucket, tag, and split
scratch. The public C boundary validates sizes, alignment, and disjointness
before calling them.

## Contracts and invariants

- Scalar bit length is public and at most 256. Scalar values are public; table
  indexing and zero-window skipping are intentionally variable-time.
- Plain paths mask to `nbits`. Decomposed paths prepare scalars modulo `r` and
  preserve exact affine identity sentinels through endomorphisms.
- Window choice, bucket count, split count, and all byte sizes use checked
  arithmetic and return zero on unrepresentable shapes.
- G2 selects four-way GLS only through 8,192 points. The public aggregate
  scratch query now applies the identical cutoff.
- Retained internal executors require nonnull explicit scratch. They allocate
  nothing and do not fall back to hidden vectors.
- Scratch belongs exclusively to one call. Concurrent calls require distinct
  storage.

## Dependency and call-path map

`blst/avx_c_api.cpp` is the production caller. It owns G1/G2/GT scratch layouts,
split placement, and validation. Direct C++ tests exercise the lower kernels
with explicit scratch. The benchmark already supplies all GT scratch. No other
production caller used the deleted allocating conveniences.

## Upstream/specification comparison

This is a CatBLST AVX implementation, not an upstream BLST file. Window
integration follows ordinary variable-time Pippenger semantics; selected
scalar decomposition and endomorphism contracts are reviewed in
`scalar_mult.h` and the accepted P4 correspondence row.

## Checklist evidence

### 1. Identity and purpose

All retained helpers have production or direct evidence callers. Dead generic
window, G2 bucket-zero, allocation, and selected GT wrappers were removed.

### 2. Contracts and invariants

Point domains, public timing, bit widths, identity sentinels, scratch shape,
ownership, and allocation behavior are recorded here and at the C declaration.

### 3. Upstream and specification lineage

Not applicable as a direct upstream file. GLV/GLS dependencies retain their
pinned VROOM and formal correspondence records.

### 4. Correctness

The audit traced window extraction, bucket indices, accumulation, infinity,
GLV/GLS/GLV4 split bounds, and scalar reduction. No arithmetic defect was
confirmed. Opposite-term cancellation now has direct G1/G2 regressions.

### 5. Memory and lifetime safety

All production paths use caller-owned live objects or byte tags. Checked size
helpers prevent wraparound. The C boundary rejects overlap before mutation.

### 6. Undefined and implementation-defined behavior

Window shifts are bounded by selected widths. Checked multiplication and range
construction reject unrepresentable counts and pointer-end wraparound.

### 7. Side channels and cryptographic hygiene

The complete file is variable-time by contract and accepts only public
scalars. It must not be reused for secret scalar multiplication.

### 8. Randomness and secret lifecycle

N/A: no randomness or secret-capable scalar contract is present.

### 9. API and ABI

Internal C++ helpers are not a supported ABI. The C surface remains unchanged;
its declarations now state disjointness, failure atomicity, allocation, and
exception behavior.

### 10. Concurrency

There is no global pool or mutable state. Scratch is exclusively caller-owned.

### 11. Resource exhaustion and denial of service

Queries fail with zero on overflow. G2 no longer requests unnecessary split
storage above its GLS cutoff. Work remains proportional to public input size.

### 12. Portability and build behavior

Portable-AVX Debug compiles and runs. Portable Release is excluded by policy;
native AVX remains the deployment profile.

### 13. Performance-sensitive structure

No bucket loop or decomposition formula changed. Overlap validation is one
fixed-count boundary check. Hidden allocations and tag fallbacks were removed.

### 14. Tests and evidence

Nineteen MSM tests and ten AVX C-contract tests pass in the focused executable.
New cases cover cutoff sizing, cancellation, overlap, failure atomicity, and
one-point aliasing.

### Specialized additions

`AUD-358` is validated by the later native Peach packet. This review neither
adds nor removes its exact-zero observations.

The implementation changes are commit `5b08edc8`; commit `4c58255b` adds only
the final GT executor contract required by the documentation ledger.

## Findings

`AUD-039`, `AUD-085`, `AUD-091`, `AUD-092`, `AUD-095`, `AUD-358`, `AUD-389`,
`AUD-390`, `AUD-391`.

## Validation

- Focused portable-Debug object compilation passed.
- Focused execution passed 37/37 tests: 19 MSM, 10 C-contract, and 8 curve.
- Structural search found no `std::vector` or removed helper in this header.
- `git diff --check` passed.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing C/C++ baseline | no retained exception path found | Does not prove MSM arithmetic |
| Clang analyzer / clang-tidy | Existing AVX triage | no mapped defect | Template and RNS bounds need direct review |
| Clippy | Existing Rust boundary review | no new binding change | C++ implementation is outside Clippy |
| ASan/UBSan | Not rerun | deferred | Focused overlap checks are direct, not sanitizer evidence |
| Fuzzing | Not run | deferred | Random MSM corpus remains milestone work |

## Residual risk and decision

Native Peach execution and before/after performance evidence close `AUD-358`.
Sanitizer coverage and fuzzing remain release gates. The file review
is complete: every retained helper and caller was traced, contracts are
explicit, and focused correctness and boundary evidence passes.
