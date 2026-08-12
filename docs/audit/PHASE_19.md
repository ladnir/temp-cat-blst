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

## Packet 19B: Retained CI and dependency automation

The reviewed repository state is
`0c7865082575e15f779debab8186bd0eca0eea7a`. Findings `AUD-347` through
`AUD-350` are validated; `AUD-351` remains open for an explicit Rust MSRV
decision.

Supported CMake CI now consumes the audited presets instead of duplicating
their flags. Linux exercises both pure C++ scalar and runtime-dispatched
assembly, Windows exercises runtime-dispatched assembly, and `macos-15`
exercises the supported native Apple Silicon preset with an explicit ARM64
check. Portable AVX remains a separate path-scoped Debug developer lane with
assembly and ADX disabled. Native AVX-512/IFMA remains a hardware release gate,
not an emulated deployment claim.

All workflows have workflow-level read-only permissions and concurrency, and
every job has a timeout. The full supported matrix runs for pull requests and
master pushes; feature-branch pushes no longer duplicate pull-request work.
Linux and Windows runner generations are fixed to Ubuntu 24.04 and Windows
2025. Apple remains fixed to macOS 15 ARM64.

Rust CI no longer restores `Cargo.lock` from cache. Registry, Git, and target
caches are keyed by OS, architecture, and the committed lockfile hash, and all
six retained Cargo test/smoke commands use `--locked`. This makes a lockfile
change invalidate the cache without permitting CI to test a cached dependency
graph in place of the commit under review.

`scripts/audit/check_ci_policy.ps1` enforces immutable external action SHAs,
fixed runner generations, permissions, concurrency, timeouts, the exact CMake
preset matrix, safe Cargo caching, locked Rust commands, Debug-only portable
AVX, and the retained weekly grouped GitHub Actions Dependabot policy. It runs
inside the lightweight audit-integrity workflow.

### Focused evidence

- checksum-verified actionlint 1.7.12 accepts all five workflows;
- all eight external action uses are immutable 40-hex commit references;
- the checkout, cache, CodeQL, and MSVC action commits resolve to their stated
  upstream tags;
- the structural CI policy gate passes with five workflows, four production
  CMake presets, and six locked Cargo commands;
- the file/upstream manifest gate passes with 810 tracked files and 63 upstream
  commits; and
- CMake parses the complete ten-preset document used by the selected matrix.

No workflow was dispatched, pushed, or otherwise executed remotely. No local
product build, benchmark, CUDA job, portable-AVX Release build, sanitizer, or
full suite ran for this policy packet.

### Residual decisions and evidence

- The crate advertises Rust 1.65, but the version-4 lockfile requires Cargo
  1.78 or newer and no MSRV lane exists. `AUD-351` requires a product decision;
  this packet does not silently raise or pretend to validate the MSRV.
- Hosted execution of the exact Linux, Windows, Apple, portable-AVX Debug,
  audit-integrity, and binding-regeneration workflows remains required.
- CodeQL remains manual and service-gated until GitHub Code Security is enabled;
  local CodeQL and analyzer evidence remains the current static-analysis gate.
- Native ADX and AVX-512/IFMA execution remains separate hardware evidence.

## Next packet after 19B

Resolve `AUD-351` before the Rust packaging review, then execute the corrected
hosted Apple/Linux/Windows CI matrix when a push is approved. The next local
file packet should audit remaining test/build support without starting a broad
test sweep.
