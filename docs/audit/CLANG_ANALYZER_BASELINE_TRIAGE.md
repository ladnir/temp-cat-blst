# Clang Static Analyzer Baseline Triage

## Scope and evidence

The first compiler-native baseline analyzes commit
`6bee5ca43f94` plus the uncommitted split-analyzer runner infrastructure. It
uses Ubuntu Clang 18.1.3 and `scan-build` with default core, Unix, C++,
dead-code, and selected security checkers plus `-analyze-headers`.

The local evidence run is
`out/audit/clang-analyzer/6bee5ca43f94-20260801T155736Z`. Separate scalar and
portable-AVX Debug profiles compiled the retained `blst_static` production
target. The profiles exclude tests, benchmarks, assembly, frozen BTx/BTE,
experimental CUDA, and removed bindings.

The scalar profile analyzed all 18 compilation/link steps. The portable-AVX
profile analyzed all 18 source compilations and produced 20 per-invocation
SARIF files; the nineteenth Ninja step was only the static-library link. The
host command wrapper timed out after twenty minutes and detached while the WSL
job continued, so the AVX log does not contain scan-build's final summary.
The helper exited and removed its guarded temporary build directory after the
remaining source analyses completed. No AVX compiler or analyzer extraction
failure occurred.

## Result disposition

Fifteen raw analyzer results reduce to five source locations and four rule
families. Repeated reports arise from the independent profiles, G1/G2 template
instantiations, and multiple retained callers.

| Checker | Raw results | Disposition |
| --- | ---: | --- |
| `core.NonNullParamChecker` | 2 | `AUD-035`. The unused `mult_window_n` points mode permits a null first point and then forms a reference through it. Current callers supply a precomputed table and pass `points == nullptr`, so the path is latent rather than currently reachable. |
| `unix.MallocSizeof` | 4 | False positive adjacent to `AUD-035`. The allocation multiplies the row width by `sizeof(ptype)` and is exactly `npoints * sizeof(*table)`, but the cast to pointer-to-array prevents the checker from recognizing the equivalence. Removing the unused allocation mode is preferable to annotating the cast. |
| `deadcode.DeadStores` | 3 | `AUD-036`. One scalar portable-backend assignment appears in both profiles; one AVX portable-intrinsic assignment appears only in the AVX profile. Both values are overwritten or never consumed. |
| `core.UndefinedBinaryOperatorResult` | 6 | False positive. Every trace ends at the same byte-select in `vec_select` after `Fp2::mul` writes both `Fp` components through fixed six-limb portable Montgomery operations. The analyzer does not preserve those complete output writes through this templated/inlined call chain and treats byte one of the output as uninitialized. Zero-initializing hot field objects solely to satisfy the analyzer would add unnecessary work. |

The ordinary Clang compilation also repeated an ignored-attribute warning for
`cc_lshift_2`: inherited `optnone` and `BLST_FORCEINLINE` attributes conflicted,
so Clang ignored the force-inline request. History shows that `optnone` was
introduced by upstream commit `a061960b` specifically to mitigate a Clang 8-11
miscompilation in `ab_approximation_n` on targets with funnel-shift
instructions. The broad `> 7` guard remained upstream, but its original commit
message identifies a bounded affected range.

Targeted Release code generation with Clang 18.1.3 showed that the broad guard
emitted `cc_lshift_2` as a separate weak, unoptimized, stack-using function in
the inverse path. Confining `optnone` to Clang 8-11 removes that symbol, retains
forced inlining on current Clang, and eliminates the warning. The affected old
compiler range now uses ordinary `inline` plus `optnone`, rather than requesting
two contradictory attributes. AppleClang 17 remains the supported CI gate
before this part of `AUD-036` is closed; no warning suppression is used.

## Proposed production changes

### AUD-035 - preserve and harden the scalar window allocation mode

`mult_window_n` currently has two CatBLST callers, both of which precompute
fixed stack tables. That current call graph does not make the points/allocation
mode disposable: reasonable source-level BLST porting is a project goal. The
user explicitly declined its removal.

Preserve both caller-supplied-table and points/allocation modes. Before changing
production code, design a compatible contract that covers `npoints == 0`, a
null first point, checked row-count multiplication, and allocation failure
without relying on an assertion that disappears under `NDEBUG`. The compact
pointer-array convention used for later contiguous points and scalars must be
retained or migrated explicitly, not accidentally simplified away.

### AUD-036 - portable-backend warning and code-generation cleanup

The review packet removes the unused low-half assignment in the portable
three-limb division tail and the unused masked `coeff1` assignment in the
portable AVX `madd52hi` emulation. Both borrow/carry computations remain
unchanged. These statements did not feed a later expression and optimizing
compilers already removed them.

The packet also limits `optnone` to the documented Clang 8-11 range and avoids
`always_inline` only within that range. Clang 18.1.3 compiled the Release scalar
library without the ignored-attribute warning, and archive comparison confirmed
that the standalone `cc_lshift_2` symbol present before the change disappeared
after it. The exact portable-plus-assembly 256- and 384-bit inverse oracles
passed 2/2. The intended direct `Madd52` equivalence test exposed the separately
tracked `AUD-037` combined-mode regression. The repaired split now keeps lane-
array `cc` emulation available beside native intrinsics while production uses
the compile-time-selected `vv` domain. The combined target compiles and
discovers 9 tests, including dedicated low/high IFMA masked-input, carry,
accumulator-wrap, and randomized cases. The Release target then built warning-
free with supported GCC 15.2 on Peach, and CTest passed all 9 direct comparisons
on AVX-512 IFMA hardware. Portable AVX field/curve integration passed 28/28. A
full Release portable-AVX rebuild was deferred after its known excessive
compilation cost manifested; no benchmark was run. AppleClang 17 CI is the
remaining validation requirement for this packet.

## CodeQL coverage relationship

The independent CodeQL run
`out/audit/codeql/6bee5ca43f94-20260801T155008Z` produced a scalar database
without extraction warnings and a portable-AVX database with the same 68
dependent-template diagnostics as the combined baseline. Clang compiled and
analyzed the affected AVX instantiations, providing compensating compiler-
native evidence without rewriting optimized templates for CodeQL.

`AUD-030` remains open because the Rust macro extraction gap and combined
feature compensation are not yet closed. Future subsystem gates must retain
focused sanitizer and scalar/AVX differential evidence for the affected group
operations; a static analyzer cannot prove arithmetic correctness or constant
time behavior.
