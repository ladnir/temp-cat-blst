# ASan and UBSan Baseline Triage

## Scope

The initial Phase 12 sanitizer baseline ran on x86-64 Ubuntu 24.04 under WSL2
with Clang/Clang++ 18.1.3, CMake 3.28.3, Ninja 1.11.1, and Debug code
generation. The source snapshot was based on commit
`290874df11a48a5a4e6c14e394d3b8a9266dd610` and included the dedicated
sanitizer target and local runner added by this packet.

The initial configuration selected portable scalar dispatch while retaining
x86 assembly providers. The remediated runner now exercises the pure C/C++
scalar implementation:

- `BLST_PORTABLE=ON`, `BLST_ASM=OFF`, `BLST_ADX=OFF`, and
  `BLST_DEFAULT_ASM=OFF`;
- native and portable AVX-512 disabled;
- tests enabled and benchmarks disabled;
- experimental CUDA, parallel MSM, and the optional VROOM vendor disabled;
- ASan and UBSan executed as separate profiles.

Sanitizer configurations build the ordinary `blst_tests` target with a
conditionally narrowed source list and without the BTx/BTE link dependency.
They therefore do not compile or execute frozen BTx/BTE tests. The sanitizer
source list covers public C and native C++
boundaries, serialization/conformance, FFT and middle product adapters,
backend selection, inverse/division arithmetic, and scalar window allocation.
This is an initial reusable baseline, not subsystem-completion evidence. The
initial evidence used the former dedicated sanitizer executable; consolidation
changes only target composition and future reproduction commands.

## Reproduction

From an x86-64 Linux checkout with Clang 18.1.3:

```sh
scripts/audit/run_sanitizers_local.sh asan
scripts/audit/run_sanitizers_local.sh ubsan
```

`all` runs both profiles sequentially. The runner snapshots tracked and
non-ignored files into Linux-native temporary storage, records tool and source
metadata beneath `out/audit/sanitizers/<run-id>/`, builds with at most the
requested worker count, and rejects a run that reports no passing tests.

The focused filter omits seven intentionally expensive stress cases:

- `BindingHeader.NativeFrBatchAffineG1DifferentialDistributions` is retained
  for its subsystem/milestone batch.
- `BindingHeader.NativeFrBatchAffineG2MatchesPippenger`, both selected-native
  MSM capacity/range cases, and both large Wbits scratch-threshold cases remain
  in their subsystem/milestone batches; individual ASan runtimes reached
  approximately 26-186 seconds in this environment.
- `InterfaceTests.Div3LimbsRandomFuzz` belongs in the bounded fuzz/regression
  corpus rather than every focused sanitizer change loop.

## Baseline result

| Profile | Run/evidence | Tests | Time | Result |
| --- | --- | ---: | ---: | --- |
| ASan + LSan | `290874df11a4-20260803T005042Z/asan/tests-focused.log` | 123 selected; 113 passed; 10 AVX-only skipped | 73.970 s | clean; no address or leak report |
| UBSan | `290874df11a4-20260803T010119Z/ubsan/tests.log` | 123 selected; 113 passed; 10 AVX-only skipped | 333.152 s | clean; no undefined-behavior report |

The ASan target built under the runner and the shortened focused command was
then executed directly from that immutable build directory. The UBSan runner
completed configuration, build, execution, and result validation end to end.
Local output under `out/` is ignored evidence and is not committed.

The ten skips are expected because the baseline selects scalar production and
does not enable AVX-512. They cover the AVX expert C API and AVX conformance
variants; their scalar counterparts execute in this target.

## AUD-041 remediation

The preferred initial configuration was pure portable C/C++ with
`BLST_ASM=OFF`. CMake exposes this option and explicitly forces
`BLST_PORTABLE=ON`, but the build fails because `blst/blst_t.hpp` directly
calls assembly-named primitives such as `from_mont_384`, `add_mod_256`, and
`mul_mont_sparse_256`, whose declarations are unavailable without the assembly
backend. A second attempt with assembly enabled and `BLST_ADX=OFF` fails at
`build/assembly.S:4`, which unconditionally rejects x86-64 without
`BLST_ADX` even though the portable assembly selection includes generic
providers.

The fix adds forced-inline `_v256` and `_v384` scalar selectors in the backend.
There is no runtime dispatch: optimized default builds call x/ADX providers,
non-ADX assembly builds call generic q providers, and no-assembly builds call
`cc_*`. Native field templates and scalar fields use this boundary;
the unconditional x86 non-ADX assembly rejection is removed. A permanent
Ubuntu Release CI matrix entry builds and tests the generic q selection.

Final remediation evidence:

| Profile | Run/evidence | Tests | Time | Result |
| --- | --- | ---: | ---: | --- |
| Pure C/C++ ASan + LSan | `e72437553704-20260803T013958Z/asan/tests-final-focused.log` | 85 selected; 75 passed; 10 AVX-only skipped | 159.026 s | clean; no address or leak report |
| Pure C/C++ UBSan | `e72437553704-20260803T021706Z/ubsan/tests.log` | 85 selected; 75 passed; 10 AVX-only skipped | 157.681 s | clean; no undefined-behavior report |
| Generic q assembly + ASan | `e72437553704-20260803T013958Z/asan-generic/tests-focused.log` | 55 selected; 55 passed | 0.564 s | clean; no address report |

The generic q archive has no unresolved x/ADX primitive and the pure-C archive
has no unresolved assembly primitive. A Clang 18 Release before/after build of
the default portable-plus-ADX configuration produced byte-identical output for
all 18 production objects and for `libblst.a` (SHA-256
`95dc8ba3ac25df1e04430b54c67d48d39dcebdc2f6ccc9921f82cb986b2e6aaa`).
No benchmark was warranted because optimized production code did not change.
A local MSVC Release library build was attempted but produced no diagnostic
before the local timeout; hosted Windows remains ordinary CI evidence rather
than a closure dependency for this configuration-specific Linux repair.

## Coverage limits

- Generated assembly is not sanitizer-instrumented; the pure-C profiles cover
  the corresponding C/C++ arithmetic paths, and the generic q profile covers
  assembly selection and integration.
- Native AVX-512 and portable AVX-512 are outside this scalar baseline.
  Portable AVX Release is prohibited in routine audit work; portable AVX Debug
  remains available for targeted subsystem work.
- Rust FFI is not exercised by this C++ executable and needs a dedicated
  Rust/FFI sanitizer packet.
- Frozen BTx/BTE, experimental CUDA, non-Apple ARM64, CHERI, and benchmarks are
  excluded by policy.
- The seven expensive tests above remain required at the appropriate subsystem
  or milestone gate.
- A clean initial baseline does not close any production file's one-by-one
  audit; sanitizer coverage expands alongside Phases 13-19.
