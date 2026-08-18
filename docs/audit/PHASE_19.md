# Phase 19: Support Files and Release Evidence

## State

Phase 19 is in progress. Packets 19A through 19F review the build, workflow,
hosted execution, and supported release evidence. Packet 19G implements the
API-complete fuzz foundation. Long fuzz discovery and the remaining packaging
and release-candidate gates remain queued.

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
optimization. Final shared installation, package configuration, and the
selectable public prefix remain Phase 20 work. `AUD-016` separately validates
the current deterministic compatibility export allowlist.

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
remained at that point. The later Peach packet closes it.

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
left `AUD-358` pending its native preprocessing performance check; the later
Peach packet completes that evidence. Generated assembly, analyzer, sanitizer,
fuzz, negative cross-toolchain,
combined-provider, packaging, and Phase 20 ABI-prefix evidence open. No product
build, benchmark, CUDA, sanitizer, fuzz campaign, or BTx/BTE work ran.

## Packet 19F: Product-owner execution closure

The current `blst/CMakeLists.txt` source predates the consolidated Peach native
AVX runs and supported Apple run `31767696988`. Peach Clang 21.1.8 and GCC 15.2
each build the native AVX-512/IFMA Release target and pass 79/79 focused tests.
Apple arm64 configures and builds the scalar-assembly target and passes 465/465
CTest cases. `AUD-342` is validated and the product target owner is complete.
Experimental CUDA, parallel MSM, non-Apple ARM64, CHERI, and MSVC ARM64 remain
outside the production claim.

## Packet 19G: API-complete fuzz foundation

The controlling design is recorded at
`b8874839ba3fc9f98891695bf509c06137446cb2`. Phase 0 is implemented through
`8c92bfdd0190ac383a28ef24ceed855eafd853bd`. Its generated ledger classifies
1,326 supported API rows: 342 C declarations, 484 reviewed native C++ methods
or free functions, 342 raw Rust FFI mappings, and 158 safe Rust functions or
methods. Private
domain carriers, reviewed refinements, typed call frames, and structural gates
prevent a campaign from silently violating a narrow arithmetic contract.

The scalar campaigns land in `6c59ed148559e43e0ac2b2b5cea798a90cac10e2`,
`68bd058fc7061eef5b49bf9a4b95671ef201bb96`, and
`83dd5b7d7c0a464cb7d17b4e9279e9229eeae854`. They cover byte boundaries,
constructed field and group arithmetic, protocols, pairing state, uniqueness,
FFT and middle product, and MSM. All nine retained corpora replay under
RelWithDebInfo ASan/UBSan and Debug assertions. Independent modulus-oracle,
domain-factory, call-frame, and seed-decoding regressions pass.

Commit `ad472e4639e8f61691a70aded9e61e6c5962bd30` adds the bounded AVX profile
and two safe Rust targets. Portable AVX Debug replay passes its retained
corpus. Pinned `cargo-fuzz` 0.13.2 replays both safe Rust targets under ASan
and instruments their native C/C++ dependency under a separate UBSan replay.
Native AVX remains manual.

Commit `486cf3bf81d6dab958fc9d29130d92a9bd3bc6b8` adds sequential runners,
resource limits, approval gates, deterministic CI replay, and manual workflow
modes. Commit `58f078b871c6a64b0f06d3574af5c57ce467d3da` emits independent LLVM
profiles and reports for all nine scalar binaries, which avoids invalid
cross-binary profile merging. Structural, workflow-policy, shell-syntax, and
workflow-lint checks pass.

Commit `9cd198a1136e765e5c3e3b5d70d2e323f015d716` closes the foundation
audit gaps. Clang now inventories native free functions as well as methods.
The structural gate binds every direct C row, unique C++ row, and safe Rust
row to executable evidence. Named success, rejection, boundary, aliasing, and
round-trip paths all execute under LLVM coverage. Manual discovery uses a
working corpus under `out/`, so automation cannot modify checked-in seeds.
The follow-up passes 14/14 scalar Debug tests, 10/10 focused
RelWithDebInfo tests, and the portable-AVX Debug replay. Rust ASan replay
reports 312 decode edges and 931 state edges; the native-UBSan replay reports
105 and 485 edges. These counts are descriptive and profile-specific.

Three observed mismatches were harness defects: MSM tile interpretation,
Debug replay time limits, and double delivery of encoded corpus seeds. Each
was corrected before accepting the replay evidence. Fuzzing found no confirmed
library defect in this packet, so no new `AUD-*` entry was created.

At source commit `23e8b2cc`, Peach Clang 21.1.8 builds the native
AVX-512/IFMA fuzz target under RelWithDebInfo ASan and UBSan. The bounded
retained-corpus replay passes 1/1 in 0.85 seconds. No long discovery, Apple
replay, benchmark, CUDA, BTx/BTE, or portable-AVX Release build ran in that
packet.

That scalar pilot also runs on Peach at `23e8b2cc` with Clang 21.1.8,
RelWithDebInfo, ASan, and UBSan. `decode` completes 152,731 executions with
1,531 edges and 2,781 features. `bytes` completes 295,566 executions with 741
edges and 1,859 features. Each target runs for 61 seconds, sequentially, with
one fuzz worker. Both artifact directories contain only `discovery.log`; no
crash, timeout, sanitizer finding, or library mismatch occurred. Generated
corpora remain under Peach's `out/` tree and do not modify checked-in seeds.

Clang 21's greedy set-cover merge initially reduces the full feature corpus to
216 additions. That packet is still too large for routine review. Repeating
the merge without counter or value-profile features produces an
edge-preserving packet of 44 decode seeds and 26 byte-conversion seeds. The
retained decode corpus preserves all 1,531 edges and 2,358 of 2,781 features;
the retained byte corpus preserves all 741 edges and 1,336 of 1,859 features.
The promoted additions total 2,012 bytes, obey both target input bounds, use
content-derived names, and replay successfully under the same ASan/UBSan
profile. A focused local sanitized replay of the promoted checked-in corpora
also passes 2/2 without rebuilding product sources.

The promotion audit found that Git's default text detection normalized five
content-addressed seeds during staging. The corpus is binary input, even when
a seed happens to resemble text. A narrow `.gitattributes` rule now disables
text transforms under `fuzz/corpus`; all 70 committed blobs match their SHA-1
filenames after restaging. This was a corpus-packaging defect, not a CatBLST
finding.

The next sequential Peach pilot covers `fields` and `groups` at source commit
`23e8b2cc` under Clang 21.1.8 RelWithDebInfo, ASan, and UBSan. `fields`
completes 1,143 executions with 2,552 edges and 5,774 features. `groups`
completes 1,595 executions with 1,988 edges and 3,825 features. Their lower
execution rates reflect the field and direct-multiplication oracles. Neither
target crashes, times out, reports a sanitizer finding, or produces an
artifact beyond `discovery.log`.

Edge-preserving set-cover minimization promotes 19 seeds for each target. The
retained field corpus preserves all 2,552 edges and 4,398 features. The
retained group corpus preserves all 1,988 edges and 2,814 features. The 38
additions total 1,037 bytes, obey both input bounds, use content-derived names,
and replay successfully under the Peach sanitizer profile. Focused local
sanitized replay also passes 2/2 without rebuilding product sources.

The first complete scalar discovery campaign runs all nine targets on Peach
at source commit `23e8b2cc`. Clang 21.1.8 uses RelWithDebInfo, ASan, and UBSan.
Fifteen fuzz workers share 16 logical CPUs, and every worker runs for 15
minutes. Each target owns a distinct corpus and artifact directory. The
controller exits successfully, and no crash, timeout, sanitizer finding,
assertion, or library mismatch produces an artifact.

The initial parallel controller placed libFuzzer job-log files from different
targets in one working directory. Those temporary logs collided. The collision
did not mix target corpora or crash artifacts, but it invalidated exact
discovery counters. This harness defect does not create an `AUD-*` finding.
The retained runner now assigns each target a private working directory.

A subsequent one-core replay executes all 2,809 discovered inputs with mutation
disabled. Every target passes under the same sanitizer profile. Edge-preserving
set-cover minimization then uses the checked-in corpus at `94aea12a` as a
protected baseline. The final corpus preserves every replayed edge:

| Target | Working inputs | Retained inputs | Edges | Retained features |
| --- | ---: | ---: | ---: | ---: |
| Decode | 690 | 75 | 1,538 | 2,545 |
| Bytes | 405 | 42 | 741 | 1,380 |
| Fields | 405 | 40 | 2,556 | 4,578 |
| Groups | 522 | 36 | 1,990 | 2,920 |
| Protocol | 83 | 21 | 2,851 | 2,994 |
| Pairing state | 157 | 8 | 3,059 | 3,665 |
| Uniqueness | 201 | 20 | 171 | 499 |
| FFT/middle-product | 97 | 20 | 3,258 | 9,300 |
| MSM | 249 | 18 | 2,372 | 5,173 |

The promotion adds 145 content-addressed seeds. The complete scalar corpus now
contains 280 inputs and 6,660 bytes. Every hash-derived filename matches its
input bytes. The standard local ASan/UBSan replay passes all 14 fuzz and
infrastructure tests in 40.08 seconds. These results provide executed boundary
and differential evidence; they do not prove arithmetic or constant-time
properties.

The first safe Rust discovery pilot runs both targets on Peach at source
commit `8f75b465`. The pinned rustc nightly 1.94.0 toolchain uses LLVM 21.1.8,
and `cargo-fuzz` is fixed at 0.13.2. GCC 15.2 builds the portable native
dependency. Each target receives eight isolated workers for five minutes per
worker. Decode completes 98,664 executions, and state completes 13,657
executions. No worker reports a crash, timeout, assertion, sanitizer finding,
or oracle mismatch.

The working decode corpus contains 854 inputs and 9,800 bytes. The state corpus
contains 80 inputs and 632 bytes. Mutation-disabled ASan replay passes 855 and
81 runs, respectively. A separate native-dependency UBSan replay passes the
same 934 inputs. Every discovered hash-derived filename matches its input
bytes, and no non-log artifact exists.

Edge-only set-cover minimization protects the exact four committed seed blobs.
The input union also includes the pilot's Windows CRLF variants. Decode retains
all 671 replayed ASan edges in 91 inputs and 609 bytes. State retains all 941
replayed ASan edges in 8 inputs and 67 bytes. The promotion adds 89 decode
seeds and 6 state seeds, totaling 613 bytes. Every addition is content-addressed
and respects its target's input bound. Standard local replay passes 101 ASan
runs and 101 native-dependency UBSan runs.

The promotion review finds that the Rust corpus lacked the binary attribute
already applied to the C/C++ corpus. Windows therefore checked out the four
retained Rust seeds with CRLF, while Git stored LF blobs. The narrow
`/bindings/rust/fuzz/corpus/** binary` rule now prevents transforms, textual
diffs, and merges. The promoted packet protects the exact Git blobs and retains
any distinct pilot behavior through content-addressed additions. This was a
corpus-packaging defect, not a CatBLST defect, so it creates no `AUD-*` item.

A local coverage-path gap audit adds one fuzz-only marker for each scalar
target that previously lacked a named path: fields, protocol, pairing state,
uniqueness, and MSM. The coverage configuration and structural self-tests
pass. A fresh local Clang 18.1.3 RelWithDebInfo ASan/UBSan replay reaches all
nine named paths with the committed 280-input corpus. The five new path counts
are 41 for fields, 22 for protocol, 9 for pairing state, 19 for uniqueness, and
19 for MSM. No crash, sanitizer diagnostic, assertion, or oracle mismatch
occurs. No product source or hot path changes.

The next packet concentrates discovery on high-level BLS interfaces. It adds a
dedicated C target for both public-key orientations. Each input drives signing,
direct verification, two-party aggregation, pairing-context merge, compressed
round trips, and altered-message rejection. A separate safe Rust target drives
the analogous `min_pk` and `min_sig` interfaces. It also covers incremental,
serialized, fast, pre-aggregated, and weighted aggregation. The existing
low-level targets remain available for deterministic replay. No fuzz-only code
enters a production path.

The Peach source starts from `153ec002`. The initial byte-exact overlay has
SHA-256 `4cb2da083e58d9d16dc3dc000313446fcfe5b6325805409663f6df82d90d7a6b`.
After oracle repair, the safe Rust target has SHA-256
`d50835cd71a6c89e4b4ad4aeae116f283a50445e5f44e3daacbaf977d41774fe`.
Clang 21.1.8 builds the C targets under RelWithDebInfo ASan and UBSan. Rustc
nightly 1.94.0 and `cargo-fuzz` 0.13.2 build the safe target under ASan. Every
build uses at most two workers.

Each target receives eight isolated workers for five minutes. Target groups run
sequentially on Peach's 16 logical CPUs. The C high-level target completes
3,725 executions, and pairing-state completes 7,471 executions. Their working
corpora contain 138 and 766 inputs. A single-process replay passes all 904
inputs, reaches 2,873 and 3,059 edges, and produces no artifact.

The first safe Rust attempt exposes a harness-oracle defect. The fast aggregate
negative check signs `message_a` but alters `message_b`. Four three-byte inputs
make the altered value equal `message_a`, so successful verification is the
correct result. The repaired target derives separate altered messages. The four
artifacts and the initial 86-input working corpus then pass sequential ASan
replay. No CatBLST source changes, and no `AUD-*` item is created.

The corrected safe Rust campaign completes 12,188 executions across eight
five-minute workers. It produces no artifact. Its 90-input combined corpus
passes a single-process replay and reaches 447 edges. Edge-only set-cover
minimization then protects the checked-in baselines and produces this retained
packet:

| Target | Working inputs | Additions | Retained inputs | Replayed edges | Addition bytes |
| --- | ---: | ---: | ---: | ---: | ---: |
| C high-level BLS | 138 | 3 | 4 | 2,873 | 26 |
| C pairing state | 766 | 5 | 13 | 3,059 | 73 |
| Safe Rust high-level BLS | 90 | 2 | 3 | 447 | 17 |

The ten promoted files total 116 bytes. Every filename matches the SHA-1 of its
bytes. The retained unions replay under the Peach sanitizer profiles. Focused
local replay also passes both C targets under ASan/UBSan and the Rust target
under ASan and native-dependency UBSan. The scalar corpus now contains 289
inputs and 6,794 bytes. The safe Rust corpus contains 102 inputs and 733 bytes.
A final local Clang 18.1.3 coverage replay reaches all ten named scalar paths.
The new high-level marker has count 5, and the pairing-state rejection marker
has count 14. The standard local C/C++ fuzz label passes 15/15 tests in 32.99
seconds. These runs provide executed high-level API evidence. They do not prove
correctness, security, or constant-time behavior.

The next high-level review extends the safe Rust target to the two-item
`verify_multiple_aggregate_signatures` flow. Both public-key orientations now
exercise successful batch verification, altered-message rejection, input-shape
rejection, invalid scalar width, and composed decode/validation failures. The
changes remain fuzz-only.

The first target revision has SHA-256
`b645a8033dd5b0a4052f5931133598646ba72b2f7933c02a99eb85421172aa7d`.
Eight isolated Peach ASan workers run it for five minutes each and complete
15,623 executions. The run produces no artifact or new corpus input, and its
90-input replay reaches 539 edges. Evidence review rejects this run as
batch-protocol closure because the target derives its weights from the signing
keys. The weights are valid and nonzero, but they do not model the API's
independence requirement.

The corrected target has SHA-256
`cd5d1d0699cb3a7205b4a1d4e8721f736d1ea062e0a6fa0c8d140988b445ced0`.
It derives verifier-weight scalars from sources that are domain-separated from
both signing keys. Focused local replay passes the three committed inputs under
ASan and native-dependency UBSan. A mutation-disabled Peach ASan replay passes
the existing 90-input corpus and reaches 555 edges.

The corrected target then runs on eight isolated Peach ASan workers for five
minutes each. The run uses Clang 21.1.8, `cargo-fuzz` 0.13.2, and at most two
build workers. Its named tmux session has a 1,800-second outer timeout. The
workers complete 15,643 executions without a crash, timeout, assertion,
sanitizer diagnostic, or oracle mismatch. They produce no artifact or new
input. A final single-process replay passes the unchanged 90-input corpus. It
reaches 555 coverage edges and 608 feature edges.

The rejected evidence concerns a fuzz harness precondition, not CatBLST, so it
creates no `AUD-*` item. The corrected run completes the planned batch-scenario
discovery budget. Neither run proves correctness, security, or constant-time
behavior.
