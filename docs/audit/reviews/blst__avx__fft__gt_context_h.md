# File Review: `blst/avx/fft/gt_context.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/avx/fft/gt_context.h` |
| Roadmap phase / queue | Phase 17 Packet 17K, after the MatrixNoK P3/P4 arithmetic owners and before final FFT milestone evidence |
| Reviewed implementation commit | `ec18911a01e097bcbdacd9f16634049b4807c998` |
| Proof/evidence rebind | `52d71c1a882fea3179fcb7b9a768026bf82f5404` |
| Reviewer | Codex audit team |
| Date | 2026-08-11 |
| Disposition | `in_scope_production` |
| Subsystem | MatrixNoK AVX GT coefficient adapter for FFT and middle product |
| Profiles | native AVX-512/IFMA, portable-AVX Debug, formal proof |
| Generation/provenance | Handwritten CatBLST adapter derived from the earlier local VROOM GT FFT integration; Apache-2.0 and Category Labs |
| Upstream references | `not_applicable: CatBLST-specific adapter around VROOM-derived MatrixNoK owners` |
| Related files | `blst/avx/gt/glv4.h`, `blst/avx/fields/matrixnok_fp12.h`, `blst/fft/ctFft.h`, `blst/fft/middle_product.h`, backend-alias and FFT tests, AVX correspondence and dispatch validators |

## Purpose and trust boundary

This header adapts selected MatrixNoK Fp12 storage to the generic FFT
coefficient policy. Its values are strictly members of the prime-order target
group GT. It provides GT identity, group multiplication, division by
cyclotomic conjugation, public Fr exponentiation, scalar/selected conversion,
and immutable precompiled public FFT-root schedules.

It must not become a general Fp12 API, validate subgroup membership, introduce
runtime backend dispatch, allocate inside butterflies, or retain process-global
mutable state. Roots, transform sizes, and GT coefficients are public. The
algorithms are intentionally variable-time in those public values.

## Contracts and invariants

- `AvxGtElement::data` is a MatrixNoK standard element representing GT. Its
  all-zero default state is only an overwrite-before-read placeholder.
- `plus` multiplies GT values, `minus` multiplies by the cyclotomic inverse,
  and `mult` exponentiates by a reduced Montgomery Fr value. None validates
  the domain in Release.
- `fr_to_limbs_le` converts reduced Montgomery Fr to a canonical four-limb
  little-endian integer. Debug checks both representation bounds.
- Prepared roots own one Fr copy and one compiled GLV4 schedule per root.
  Construction may allocate and throws normally on allocation failure.
- `prepare_inverse_scale()` requires a nonempty table and derives `1/n` from
  the owned table size; callers cannot supply a divergent scale length.
- Prepared public fields are expert-visible but immutable after construction.
  Mutation, movement, or destruction while scoped is outside the contract.
- `ScopedRoots` is synchronous, noncopyable, nonmovable, and lvalue-only. It
  may nest on one thread and restores the prior thread-local pointer during
  normal or exceptional scope exit.
- A prepared root schedule is selected only for exact pointer identity into
  the active root vector. An equal but unrelated Fr follows the stack-compiled
  fallback. The separately prepared scale is selected by public value equality.

## Dependency and call-path map

`blst/fft/middle_product.h` is the sole production owner. It constructs inverse
roots, prepares the self-derived scale, enters `ScopedRoots`, and invokes the
generic full/truncated FFT with `AvxGtCoeffCtx`. Direct FFT and backend-alias
tests cover the same paths. Public C and safe Rust middle-product calls reach
this header only through the reviewed native owner and its exception-translating
C boundary; this header exports no C symbol and throws only during precompute
allocation.

Arithmetic delegates to the selected MatrixNoK Fp12 ring and accepted GLV4
compiler/executor. Scalar/default bridges use the reviewed direct typed Fp12
conversion. Native AVX feature gating remains owned by build/dispatch policy,
not this header.

## Correctness, safety, and performance review

Every declaration and branch was reviewed. `AUD-028` and `AUD-029` already
closed temporary lifetime and unrelated-pointer ordering risks. The exact
pointer-selection sequence never subtracts unrelated pointers and checks the
derived index plus identity before access. Vector construction is exception
safe; no prepared pointer is published by the constructor. No secret, entropy,
recursion, unbounded stack object, foreign callback, or C exception boundary is
owned here.

`AUD-313` records and repairs the missing GT-only, reduced-Fr, timing, allocation,
immutability, and scoped-lifetime contracts. `AUD-314` removes the dead 32-byte
canonical limb vector and scale field: canonical limbs now exist only in fixed
stack arrays while their schedules compile. This removes one allocation and 32
bytes per root plus 32 bytes per precompute. `AUD-315` makes inverse scale a
self-derived property of the root owner.

The Release `mult` arithmetic and selection order are unchanged. New invariant
checks compile out with assertions, return-type spelling is explicit, and ring
owners remain constant-initialized. No benchmark was warranted: the hot path is
unchanged while preparation strictly removes allocation and retained storage.

## Formal and independent evidence

Accepted row `AVX-P3-007` proves the compiled GT executor for every canonical
Fr scalar below `r`, including schedule and physical bounds and the FFT scalar
corollary. `AVX-P4-002` binds selected production ownership. Commit `52d71c1`
rebinds both rows to the exact reviewed header; 21 rows and 61 mapped hashes
validate. The change does not alter `glv4.h`, any proved arithmetic expression,
or the compiled schedule consumed by `mult`.

Scalar/AVX stage replay checks group operations after every FFT stage, compares
prepared root and inverse-scale schedules with the fallback compiler, and
compares full and truncated inverse transforms. Historical Peach Clang 21.1.8
native Release evidence covers the unchanged arithmetic and owner initialization.

## Automated evidence

| Tool/profile | Result | Blind spot / residual |
| --- | --- | --- |
| CodeQL baseline | `AUD-028` and `AUD-029` validated | Does not prove GT algebra or compiler constant-time behavior |
| MatrixNoK dispatch validator | Passes; pins selected owners and forbids dead caches/helpers | Structural only |
| AVX correspondence validator | 21 accepted rows and 61 hashes pass | Compiler/ISA lowering remains P0 |
| MSVC portable-AVX Debug | Target rebuilt with two workers; four focused tests pass | Emulates AVX and is not timing evidence |
| Native AVX-512/IFMA | Historical accepted arithmetic evidence retained | Exact post-cleanup native release rerun deferred to milestone gate |
| ASan/UBSan and fuzz | Not run in this packet | Remain Phase 17/19 release gates |

## Findings

`AUD-028`, `AUD-029`, `AUD-216`, `AUD-220`, `AUD-229`, `AUD-313`, `AUD-314`,
and `AUD-315`.

## Validation

- `cmake --build out/build/x64-Debug-Scalar --target blst_tests --parallel 2`
  passes; frozen BTx/BTE rebuilt unchanged and emitted its existing nodiscard
  warning.
- Four focused tests pass: invalid prepared input death checks, prepared and
  fallback root/scale equivalence, full inverse FFT, truncated inverse-stage
  replay, and scalar/AVX middle-product agreement.
- `python scripts/audit/validate_matrixnok_dispatch.py` passes.
- `python scripts/audit/validate_avx_correspondence.py` passes all 21 rows and
  61 hashes.
- Structural search finds no retained limb vector/scale, parameterized scale
  setter, stale AVX52 test spelling, or test-only conversion helper.
- No benchmark, CUDA, portable-AVX Release, full suite, sanitizer, remote job,
  Apple run, or BTx/BTE audit/change was performed.

## Residual risk and decision

The file is complete for hand review and accepted proof binding. Native
post-cleanup compiler/ISA lowering, sanitizer, fuzz, Apple, and broader release
evidence remain milestone gates. Expert mutation and forged/scoped lifetime
misuse remain narrow C++ precondition violations. GT membership is established
by reviewed callers rather than checked here.
