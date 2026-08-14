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
| `core.NonNullParamChecker` | 2 | `AUD-035`, validated. The retained points mode now rejects a null first point before pointer formation while preserving later null entries as the compact contiguous convention. |
| `unix.MallocSizeof` | 4 | `AUD-035`, validated. The retained allocation mode now checks `npoints * sizeof(*table)` for overflow and uses that exact row-sized expression, removing both the ambiguity and the unsafe multiplication without deleting BLST-porting functionality. |
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
two contradictory attributes. Native Apple Silicon CI at commit `0b7ec47`
subsequently compiled this path with AppleClang 17 without the former
conflicting-attribute warning; no warning suppression is used.

## Production dispositions

### AUD-035 - preserve and harden the scalar window allocation mode

`mult_window_n` has two CatBLST production callers, both of which precompute
fixed stack tables. The points/allocation mode remains available for reasonable
source-level BLST porting, as the user required. The function now returns
status and leaves the output unchanged when the point/scalar preconditions,
nonzero sizes, checked row-count multiplication, or allocation fail. The
caller-populated and precomputed table modes remain intact, as do later null
entries that select contiguous points and scalars.

Allocation injection is a compile-time policy used to test failure
deterministically; it introduces no function pointer, allocator object, or
runtime dispatch. Four dedicated G1/G2 tests cover allocated,
caller-populated, and precomputed tables, compact pointers, every rejection,
overflow, and allocation failure. Together with the existing G1 GLV and G2 GLS
multiplication/signing paths, 10/10 focused tests passed under MSVC Release,
Linux ASan with leak detection, and Linux UBSan. An optimized GCC 13 object
contained the two production `mult_glsv` symbols but no `mult_window_n`,
`malloc`, or `free` symbol/reference. No benchmark was run because the
precomputed hot-path code retained no allocation or out-of-line adapter.

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
compilation cost manifested; no benchmark was run. GitHub Actions run
`30754724252` supplied the final supported-Apple gate on native ARM64:
AppleClang 17 configured and built with `BLST_ALLOW_UNTESTED_ARM64=OFF`, all
336 C/C++ tests passed, and the Rust default and full-binding invocations
passed. `AUD-036` is validated.

## CodeQL coverage relationship

The independent CodeQL run
`out/audit/codeql/6bee5ca43f94-20260801T155008Z` produced a scalar database
without extraction warnings and a portable-AVX database with the same 68
dependent-template diagnostics as the combined baseline. Clang compiled and
analyzed the affected AVX instantiations, providing compensating compiler-
native evidence without rewriting optimized templates for CodeQL.

`AUD-030` is validated by explicit compensation. The combined
`full-rust-bindings,serde-secret` Clippy/MSRV profiles compile the macro-expanded
Rust surface. Focused scalar/AVX differentials, completed AVX owner reviews,
and formal correspondence cover the affected group operations. Future gates
must retain those controls because static analysis does not prove arithmetic
correctness or constant-time behavior.
