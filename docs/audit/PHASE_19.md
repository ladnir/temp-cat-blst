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
`a71068a47c361c9a1edae4a706b75435023756f6`. Findings `AUD-347` through
`AUD-351` are validated.

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
seven retained Cargo test/smoke/check commands use `--locked`. This makes a lockfile
change invalidate the cache without permitting CI to test a cached dependency
graph in place of the commit under review.

The formerly ignored lockfile is now a tracked CI reproducibility input. Rust
1.78 is the declared and tested MSRV because it is the first Cargo release that
accepts lockfile version 4. A separate Ubuntu job installs exactly Rust 1.78.0
and checks the combined `full-rust-bindings`, `portable`, and `serde-secret`
feature surface with two workers and a 30-minute bound. Weekly grouped Cargo
Dependabot proposals must pass this lane before review or merge.

`scripts/audit/check_ci_policy.ps1` enforces immutable external action SHAs,
fixed runner generations, permissions, concurrency, timeouts, the exact CMake
preset matrix, safe Cargo caching, locked Rust commands, Debug-only portable
AVX, the exact Rust 1.78/v4-lock contract, and the retained weekly grouped
GitHub Actions and Cargo Dependabot policies. It runs
inside the lightweight audit-integrity workflow.

### Focused evidence

- checksum-verified actionlint 1.7.12 accepts all five workflows;
- all eight external action uses are immutable 40-hex commit references;
- the checkout, cache, CodeQL, and MSVC action commits resolve to their stated
  upstream tags;
- the structural CI policy gate passes with five workflows, four production
  CMake presets, and seven locked Cargo commands;
- Rust/Cargo 1.78 checks the complete supported feature combination in 29.85
  seconds locally with one build worker;
- the file/upstream manifest gate passes with 818 tracked files and 63 upstream
  commits; and
- CMake parses the complete ten-preset document used by the selected matrix.

No workflow was dispatched, pushed, or otherwise executed remotely. No local
product build, benchmark, CUDA job, portable-AVX Release build, sanitizer, or
full suite ran for this policy packet.

### Residual evidence

- Hosted execution of the exact Linux, Windows, Apple, portable-AVX Debug,
  Rust 1.78, audit-integrity, and binding-regeneration workflows remains
  required.
- CodeQL remains manual and service-gated until GitHub Code Security is enabled;
  local CodeQL and analyzer evidence remains the current static-analysis gate.
- Native ADX and AVX-512/IFMA execution remains separate hardware evidence.

## Next packet after 19B

Execute the corrected hosted Apple/Linux/Windows/Rust CI matrix when a push is
approved. The next local
file packet should audit remaining test/build support without starting a broad
test sweep.

## Packet 19C: First hosted execution

The first hosted execution of Packet 19B established substantial positive
evidence: audit integrity, binding regeneration, Rust 1.78 MSRV, current Rust
on Linux/Windows/macOS, Linux pure-C++ scalar, and Linux assembly all passed.
The supported Apple product archive also built, and the Windows C++ suite
passed 449 of 450 tests.

Three portability defects prevented complete gates and are recorded as
`AUD-352` through `AUD-354`. Commit `d69be47` repairs the invalid C++20 test
comparison signature, normalizes the generated-proof source hash across LF and
CRLF checkouts, and matches the fixed-array AVX batch template parameter to
`std::size_t`. Focused Clang, GCC, and MSVC compilation plus an LF/CRLF
generator differential pass. No product arithmetic, runtime dispatch, or hot
loop changed; the AVX repair is a compile-time type correction only.

The exact hosted revision has not yet been rerun, so Apple, Windows, and the
portable-AVX developer lane remain release-evidence gaps. No benchmark, CUDA,
BTx/BTE work, full local suite, sanitizer sweep, or broad parallel build ran.

## Packet 19D: Second hosted execution

The second hosted execution confirmed the three Packet 19C repairs: Apple
compiled the FFT regression, Windows accepted the line-ending-independent Lean
snapshot, and GCC compiled the portable fixed-batch API. Audit integrity and
all Rust jobs passed again. Linux scalar, Linux assembly, and Windows completed
their product suites with only the omitted proof-ledger hash failing.

The deeper portable-AVX Debug execution exposed `AUD-356` through `AUD-358`:
an incorrectly registered combined-provider executable, a test that passed an
intentional Miller failure marker to final exponentiation, and a real G1/G2
MSM identity-representation defect in GLV/GLS preprocessing. Apple separately
exposed the `AUD-355` timer-test expression defect.

Commit `ac60a7d` repairs these boundaries. Identity preservation uses exact
RNS sentinel checks in the variable-time public-scalar preparation path; it
does not export coordinates to scalar arithmetic, allocate, or change the
non-identity endomorphism formulas. Both an in-batch and tail identity are now
covered for G1 and G2. The combined compatibility executable is registered
only when both native and portable providers exist, malformed Miller output is
checked without violating final-exponentiation preconditions, and the exact
proof hashes are rebound.

Focused Clang/GCC/MSVC translation-unit compilation, portable-only CTest graph
inspection, proof correspondence, dispatch validation, and CI structural
policy pass. No full local build or suite, benchmark, CUDA, sanitizer, native
AVX execution, or BTx/BTE audit ran. The exact hosted rerun remains required.

The third hosted execution made Linux scalar, Linux assembly, Windows, audit
integrity, Rust MSRV, and current Rust on every supported operating system
green. Apple alone showed that explicitly constructing `std::string_view` was
still part of the `noexcept` expression and is not required by the standard to
be nonthrowing. Commit `b439a49` uses unevaluated `declval` arguments, so the
assertions now test only `TimingStats` members. Focused Clang and MSVC
compilation passes; the exact Apple rerun remains required. Portable AVX was
still executing when this correction was prepared.

Final hosted evidence is green. Supported run `31665323681` passes Apple
Silicon C++ and Rust, Linux scalar and assembly, Windows, Rust MSRV, and every
focused test. Audit-integrity run `31665323679` passes both manifest and binding
policy gates. Portable AVX Debug run `31664863148` configures, builds, and
passes its complete developer-emulation test step. This closes the hosted
execution gaps for `AUD-352` through `AUD-358`; only the separately tracked
native performance evidence for the `AUD-358` preprocessing observation
remains.

## Packet 19E: Hosted-evidence reconciliation

The green runs are reconciled into the per-file records rather than left only
as phase prose. The supported and developer CI workflow owners,
`blst/common/cpuid.cpp`, and `blst/common/bytes.h` move from `in_review` to
`complete`. `AUD-057`, `AUD-058`, `AUD-062`, `AUD-063`, `AUD-066`, and
`AUD-352` through `AUD-357` move to `validated`. The assembly declaration
owner, root/build CMake owners, presets, Fp2, retained FFT, and two affected
test owners remain `in_review`, but their stale Apple or portable-execution
gaps are removed.

This packet does not claim native AVX execution from portable emulation. It
leaves `AUD-358` fixed pending its native preprocessing performance check and
leaves generated assembly, analyzer, sanitizer, fuzz, negative cross-toolchain,
combined-provider, packaging, and Phase 20 ABI-prefix evidence open. No product
build, benchmark, CUDA, sanitizer, fuzz campaign, or BTx/BTE work ran.
