# CodeQL Baseline Triage

## Scope and evidence

This record triages the first pinned CodeQL 2.26.2 baseline at commit
`e026202ed606ce29ca94f5ceefa55913ebb249de`. C/C++ used
`codeql/cpp-queries@1.8.0:codeql-suites/cpp-security-and-quality.qls` and Rust
used
`codeql/rust-queries@0.1.39:codeql-suites/rust-security-and-quality.qls`.

The C/C++ run is `e026202ed606-20260801T005352Z`. The clean Rust run is
`e026202ed606-20260801T011000Z`. Raw SARIF and databases remain ignored local
evidence; durable dispositions are recorded here and in `FINDINGS.md`.
Benchmarks, frozen BTx/BTE, and experimental CUDA were not analyzed.

## C/C++ result disposition

All 958 emitted results are accounted for. Eight results map to five findings;
the other 950 are documented non-defects or style observations.

| CodeQL rule | Count | Disposition |
| --- | ---: | --- |
| `cpp/commented-out-code` | 921 | `AUD-038`. Removed the 171 abandoned implementations, stale debug fragments, and code-shaped interface comments across 30 retained production files. The 750 reciprocal/square-root addition-chain generation traces remain as intentional arithmetic provenance. |
| `cpp/unused-local-variable` | 20 | False positives. All are compile-time non-type template arguments or type witnesses used through template arguments or `decltype`; they create no dead runtime work. |
| `cpp/large-parameter` | 4 | `AUD-027`. One hot 768-byte AVX point parameter is actionable. Three large `IntRNS2System` constructor parameters are initialization-only efficiency cleanups. |
| `cpp/unused-static-function` | 4 | False positives caused by combining scalar and portable-AVX builds in one database. Every function has direct retained call sites; no helper should be removed. |
| `cpp/sizeof-side-effect` | 2 | Non-defects. `sizeof` does not evaluate `generator()`. Prefer `sizeof(out.Z)` for clarity with no code-generation change. |
| `cpp/stack-address-escape` | 1 | Current call sites are lifetime-safe scoped RAII. `AUD-028` adds an rvalue guardrail for future callers without copying or allocating prepared roots. |
| `cpp/incomplete-parity-check` | 1 | Non-defect. The loop variable starts at one and only increases. `i % 2 != 0` is an optional equivalent cleanup. |
| `cpp/poorly-documented-function` | 1 | The documentation score is not a correctness result. Reviewing the function exposed the independent public-boundary defect `AUD-025`. |
| `cpp/suspicious-sizeof` | 1 | `AUD-024`: confirmed pointer-sized short copy into an affine signature. |
| `cpp/rule-of-two` | 1 | Non-defect. `Factor` retains an implicitly generated copy-assignment operator; the explicit default copy constructor is redundant and may be removed under the rule of zero. |
| `cpp/array-in-interface` | 1 | Non-defect internal lambda. Its array parameters decay to pointers, but only two local four-limb arrays are passed. A fixed-array reference is optional style cleanup, subject to CUDA compilation compatibility. |
| `cpp/missing-header-guard` | 1 | `AUD-026`: confirmed repeat-inclusion defect in the scalar MSM header; the directly included bulk-addition header has the same issue. |

### Exact grouped false-positive evidence

- `map_to_G.h:709,729`: `nelem` is the non-type template argument to
  `hash_to_field<nelem>`.
- The retained AVX unused-local diagnostics in `scalar_mult.h`,
  `change_base.h`, `bounded_element.h`, `fp2.h`, and `ring_element.h` carry
  compile-time bounds, sizes, or deduced types. The former `bounded_ring.h` and
  `elementwise_reduce.h` owners were later removed as unreachable.
- `middle_product_c_api.cpp:50,157,171` is called at lines 240/266, 234, and
  260. `sign/aggregate.cpp:37` is called at line 605. The merged database
  retained mutually exclusive scalar and AVX profile variants.
- `point.cpp:597,1292` uses an unevaluated `sizeof` operand; no generator call
  or side effect occurs.
- `avx/curve/scalar_mult.h:36` operates only on a positive increasing loop
  index.

## Confirmed and adjacent findings

### AUD-024 - affine signature short copy

`Point<Field>::sign2` passes `sizeof(sig)` to `vec_copy`, so a 64-bit build
copies eight bytes rather than the complete `AffinePoint<Field>`. The sibling
`to_pk2` implementation correctly uses `sizeof(*PK)`. G1 and G2 callers can
receive an invalid, mostly stale affine signature.

Resolution: `sizeof(*sig)` now copies the complete affine object without
changing the existing prefix-layout operation. Direct G1/G2 regressions use
nonzero-prefilled destinations and cover affine-only, serialized-only, and
combined outputs. Focused scalar and Linux ASan runs passed.

### AUD-025 - FFT suffix shape and scratch arithmetic

`interpolate_fft_domain_suffix_with_scratch` does not require `n <= fft_len`.
It writes `domain_work[coeff - 1]` for `coeff` through `n`, although
`domain_work` contains only `fft_len` elements. Later quotient-index arithmetic
also relies on the same missing invariant. The public C API can therefore
perform an out-of-bounds write for `n > fft_len`.

The scratch query independently computes `3 * fft_len + 3 * n + 2` without
checked arithmetic. A wrapped result can allocate or accept undersized scratch
and then partition spans beyond it.

Resolution: the C and safe Rust boundaries reject `n > fft_len`; the scratch
query and adjacent polynomial growth use checked arithmetic; and invalid
shapes fail before output writes. Focused C++ boundary tests passed under Linux
ASan, while the default and portable Rust FFT module suites passed 11/11 with
suffix-oracle, oversubscribed-domain, and short-scratch coverage. The same
review also produced and closed `AUD-031` and `AUD-032` for prepared-domain
authentication and FFT size growth.

### AUD-026 - repeat-includable scalar MSM headers

`blst/scalar/curve/multi_scalar.h` defines templates without a header guard or
`#pragma once`. Its direct dependency `bulk_addition.h` has the same defect.
Current translation units include them once, but a supported consumer can
trigger redefinitions through repeated or mixed direct/wrapper inclusion.

Resolution: both headers now use `#pragma once`. A focused scalar translation
unit includes each directly more than once and through
`blst/curve/multi_scalar.h`; the dedicated scalar build and focused test pass.

### AUD-027 - avoidable AVX object copies

`GLSMult::minus_psi` takes a 768-byte projective point by value although it
does not modify it. The default 16-entry table invokes the function 48 times
per multiplication across three psi transforms. Passing by const reference
removes ABI-level copies when inlining or copy elimination does not.

The three construction-only `IntRNS2System` parameters are 2,048, 704, and 704
bytes. Const references can avoid an extra copy for lvalue callers while the
object retains its owned member copy. This is initialization-path cleanup, not
a runtime-indirection redesign.

Resolution: the hot point and three construction-only inputs are const
references while all members remain owned values. The portable AVX-selected
build and G2 Pippenger/scalar-multiplication oracle pass. No benchmark is
required for this signature-only removal of copies.

### AUD-028 - scoped prepared-root lifetime guardrail

All current `AvxGtCoeffCtx::ScopedRoots` call sites are safe: the prepared
roots lvalue outlives the guard, use is synchronous, exception unwinding
restores the prior pointer, nesting is preserved, and storage is thread-local.
Copying or allocating the prepared vectors would add unacceptable FFT setup
overhead.

The const-reference constructor could nevertheless bind a temporary that dies
at the end of the construction full-expression. Resolution: prepared-root
rvalue construction and guard copying/moving are deleted, while lvalue
construction and destruction are explicitly `noexcept`. Compile-time ownership
assertions pass. A focused nested-scope exception test proves restoration of
the outer prepared table and final fallback behavior.

### AUD-029 - unrelated pointer ordering in AVX root lookup

`AvxGtCoeffCtx::mult` uses built-in relational comparisons between `&b` and
the prepared vector's pointers. Ordering unrelated object pointers is not a
portable membership test and can admit an invalid subtraction under a
non-address-order implementation.

Resolution: lookup uses `std::less<const Fr*>`, constant-time integer
offset/alignment checks, and exact pointer identity before indexing. It never
subtracts potentially unrelated pointers and adds no allocation or scan.
Prepared-table selection, unrelated same-value fallback, nested restoration,
and scalar differential checks pass.

## Extraction coverage

`AUD-030` controls the analyzer coverage gap. C/C++ analysis emitted 68
partial extraction diagnostics: eight root errors at structured bindings of
dependent variadic `std::array` results and 60 cascades. Sixty-seven point into
`blst/avx/curve/ec.h` and one into `blst/avx/rns/bounded_ring.h`. The real GCC
build succeeds, both headers are otherwise extracted, and database telemetry
reports about 99.91% of lines and 99.98% of functions without extractor
errors. The blind spot is specifically the instantiated AVX G1/G2 add,
mixed-add, double, and dadd paths.

The clean Rust run emitted zero alerts and scanned 14 of 14 source files. It
reported six equivalent `format!` macro-expansion warnings in `lib.rs` for
three serde error-formatting sites instantiated by both signature variants.
They contain no dynamic format string or cryptographic control flow, but the
file remains imperfectly extracted.

Do not rewrite optimized AVX templates solely to accommodate an extractor.
The next analyzer packet should use separate scalar and portable-AVX
databases, retain explicit profile labels, and cover the AVX group operations
with compiler-native Clang analysis, focused sanitizers, and scalar
differential tests. File a minimized extractor reproducer and recheck it when
CodeQL is upgraded. Rust compensation must include combined
`full-rust-bindings serde-secret` compile/tests plus min-pk/min-sig serde
round-trip and malformed-input cases; separate feature jobs do not instantiate
that combination.

## Proposed implementation order

1. `AUD-024` and `AUD-025`: correctness and memory-safety fixes with focused
   regression/sanitizer checks.
2. `AUD-026`, `AUD-028`, and `AUD-029`: header and lifetime/portability
   hardening. Completed.
3. `AUD-027` and no-codegen style cleanups: AVX performance/API cleanup.
   Completed for the reported large parameters.
4. `AUD-030`: split analyzer profiles and add compensating Clang/Rust evidence.

The user approved and validated the `AUD-024`/`AUD-025` correctness packet and
the `AUD-026` through `AUD-029` focused hardening packet. The later compiler-
native, Rust feature-matrix, focused runtime, and proof evidence validates
`AUD-030` without rewriting optimized templates for the extractor.

## Split current-tree confirmation

The current-tree run `6bee5ca43f94-20260801T155008Z` creates independent
scalar and portable-AVX C/C++ databases:

- scalar emitted 906 results and no extraction warnings;
- portable AVX emitted 950 results and the same 68 dependent-template
  extraction warnings as the combined baseline; and
- the prior large-parameter and missing-header-guard alerts are gone.

The remaining non-comment CodeQL results are the already documented
compile-time/template false positives and style observations, the intentional
scoped prepared-root pointer controlled by the hardened `AUD-028` guard, and
the positive-loop parity expression. Two platform-sensitive FFT overflow
guards are redundant on the analyzed 64-bit target but remain required for
supported-width-independent sizing logic.

Compiler-native compensation is recorded in
`CLANG_ANALYZER_BASELINE_TRIAGE.md`. It successfully traversed the AVX
instantiations that CodeQL could not extract and opened `AUD-035` and
`AUD-036` for separate production review.

## Commented-code cleanup confirmation

`AUD-038` resolves the actionable portion of the CodeQL commented-code batch.
Disabled alternate implementations, stale debug output, duplicate wrappers,
and code-shaped interface descriptions were removed or rewritten across 30
retained production files. The wholly commented-out and unused
`test/BasicFp.h` was deleted separately; benchmarks, frozen BTx/BTE, and
experimental CUDA were not included in this pass.

The current-tree run `01eee6b03950-20260802T223754Z` completed both split
C/C++ profiles. Each profile now reports exactly 750
`cpp/commented-out-code` results:

- 377 in generated `blst/scalar/fields/recip-addchain.h`; and
- 373 in generated `blst/scalar/fields/sqrt-addchain.h`.

No result remains in any other file. These comments record the generated
addition-chain schedule and are retained as arithmetic provenance until a
repository-owned generator and reproducibility gate can replace that value.
This is an explicit, narrow disposition rather than a global suppression.
