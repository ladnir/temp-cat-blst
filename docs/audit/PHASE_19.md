# Phase 19: Support Files and Release Evidence

## State

Phase 19 is in progress. Its first packet reviews the CMake target graph,
configuration contracts, presets, and test dependency boundary. Tests,
workflows, vectors, remaining repository policy, packaging, analyzers,
sanitizers, fuzzing, and the supported release matrix remain queued.

BTx/BTE remains frozen, benchmarks remain excluded from file-by-file audit,
and CUDA remains experimental. This packet did not edit or audit those source
families and ran no benchmark, CUDA build, full suite, remote job, or broad
parallel compilation.

## Packet 19A: CMake and preset contracts

The reviewed implementation is `20c075a5000b8f7c29bb06a714112b1c8f0f9742`.
Findings `AUD-340` through `AUD-346` cover:

- implicit host-native code generation and public compile-option leakage;
- misleading or inert options, benchmark/test defaults, and contradictory
  backend selections;
- unsupported architecture and native-ISA admission;
- ordinary C++20 versus experimental CUDA/parallel-MSM C++23 requirements;
- mutable GoogleTest fetching;
- stale/developer-only preset exposure; and
- CUDA configuration pollution while CUDA was disabled.

Ordinary scalar/assembly builds no longer use `-march=native` or
`-mtune=native`. ADX/BMI2 and AVX-512/IFMA are explicit deployment contracts.
The selected AVX native C++ facade propagates its exact C++ ISA requirement;
ordinary scalar consumers receive only the required C++20 language feature.

Unavailable shared-library and removed VROOM-vendor modes fail configuration
instead of being silently ignored. `BLST_FAST_BUILD` remains the historical
CUDA-only choice to reduce compilation time at the cost of runtime
optimization. Final shared installation, exports, package configuration, and
the selectable public prefix remain coupled to `AUD-016` and Phase 20.

## Focused evidence

- CMake parses ten visible presets: Linux scalar/assembly/native/AVX,
  Windows scalar/assembly/native, and Apple Silicon. No visible preset enables
  benchmarks, CUDA, or portable AVX emulation.
- A fresh MSVC 19.50 scalar Release configuration emits C++20, no C++23/latest
  mode, no `/arch`, and no host-native tuning for `blst_static`.
- A serial scalar `blst_static` build and linked downstream
  `add_subdirectory` consumer pass. The consumer inherits C++20 but no private
  warning, `-fwrapv`, optimization, or CPU flags; embedded tests and benchmarks
  are both off.
- A fresh runtime-dispatched x86-64 assembly Release archive builds serially,
  including the checked-in Win64 assembly and capability owner.
- Shared mode, CUDA-fast without CUDA, CUDA-MSM without CUDA, removed VROOM
  vendor mode, default assembly without assembly, and native AVX under the
  unsupported local MSVC compiler all fail configuration.
- A clean test configuration verifies the GoogleTest 1.14 tag ZIP against
  SHA-256 `1f357c27ca988c3f7c6b4bf68a9395005ac6761f034046e9dde0896e3aba00e4`.
- The experimental parallel-MSM configuration requests C++23/latest; ordinary
  production remains C++20.

## Residual evidence

The following are required before this packet can contribute to a release
claim:

- native Linux scalar/assembly and AVX-512/IFMA configuration and execution;
- hosted Apple Silicon preset execution;
- explicit cross-toolchain rejection evidence for 32-bit, unknown, universal
  Apple, unqualified ARM64, and CHERI configurations;
- an experimental CUDA configure/build pass owned by the CUDA scope, not the
  production gate;
- installed C/C++ consumers, shared-library/export design, and Windows symbol
  allowlist closure; and
- the Phase 19 CI, analyzer, sanitizer, fuzz, coverage, package, and complete
  supported-matrix gates.

## Next packet

Audit GitHub Actions and dependency automation against the corrected CMake
contracts before executing hosted release evidence. The workflow packet must
keep portable AVX developer coverage separate from production deployment
claims and must not add CUDA to the supported release matrix.
