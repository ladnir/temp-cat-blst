# ASan and UBSan Baseline Triage

## Scope

The initial Phase 12 sanitizer baseline ran on x86-64 Ubuntu 24.04 under WSL2
with Clang/Clang++ 18.1.3, CMake 3.28.3, Ninja 1.11.1, and Debug code
generation. The source snapshot was based on commit
`290874df11a48a5a4e6c14e394d3b8a9266dd610` and included the dedicated
sanitizer target and local runner added by this packet.

The configuration deliberately selects the portable scalar dispatch while
retaining the x86 assembly providers required by the current build:

- `BLST_PORTABLE=ON`, `BLST_ASM=ON`, `BLST_ADX=ON`, and
  `BLST_DEFAULT_ASM=OFF`;
- native and portable AVX-512 disabled;
- tests enabled and benchmarks disabled;
- experimental CUDA, parallel MSM, and the optional VROOM vendor disabled;
- ASan and UBSan executed as separate profiles.

`blst_sanitizer_boundary_tests` is a small sanitizer-only executable. It does
not compile the ordinary all-tests executable and therefore does not compile
or execute frozen BTx/BTE tests. Its source list covers public C and native C++
boundaries, serialization/conformance, FFT and middle product adapters,
backend selection, inverse/division arithmetic, and scalar window allocation.
This is an initial reusable baseline, not subsystem-completion evidence.

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

The focused filter omits two intentionally expensive cases:

- `BindingHeader.NativeFrBatchAffineG1DifferentialDistributions` is retained
  for its subsystem/milestone batch; it took approximately 95 seconds alone
  under ASan in the baseline environment.
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

## AUD-041 configuration finding

The preferred first configuration was pure portable C/C++ with
`BLST_ASM=OFF`. CMake exposes this option and explicitly forces
`BLST_PORTABLE=ON`, but the build fails because `blst/blst_t.hpp` directly
calls assembly-named primitives such as `from_mont_384`, `add_mod_256`, and
`mul_mont_sparse_256`, whose declarations are unavailable without the assembly
backend. A second attempt with assembly enabled and `BLST_ADX=OFF` fails at
`build/assembly.S:4`, which unconditionally rejects x86-64 without
`BLST_ADX` even though the portable assembly selection includes generic
providers.

This is recorded as proposed low-severity finding `AUD-041`. The sanitizer
packet does not change production dispatch or assembly. Until the finding is
resolved, sanitizer instrumentation reaches the C/C++ code around arithmetic
but cannot instrument the generated assembly primitives themselves.

## Coverage limits

- Generated assembly is not sanitizer-instrumented.
- Native AVX-512 and portable AVX-512 are outside this scalar baseline.
  Portable AVX Release is prohibited in routine audit work; portable AVX Debug
  remains available for targeted subsystem work.
- Rust FFI is not exercised by this C++ executable and needs a dedicated
  Rust/FFI sanitizer packet.
- Frozen BTx/BTE, experimental CUDA, non-Apple ARM64, CHERI, and benchmarks are
  excluded by policy.
- The two expensive tests above remain required at the appropriate subsystem
  or milestone gate.
- A clean initial baseline does not close any production file's one-by-one
  audit; sanitizer coverage expands alongside Phases 13-19.
