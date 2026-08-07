# File Review: `blst/c_api.cpp`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/c_api.cpp` |
| Roadmap phase / queue | Phase 18 public C boundary, Wave 1 |
| Reviewed commit | `02f5e5ff078b251ca3b5f86f4f83ea64af664bc7` (implementation); public declaration companion `aa4c7ed187f5a63496f0e94f9284108cbd46c7f1` |
| Reviewer | Codex audit team |
| Date | 2026-08-06 |
| Disposition | in-scope production, complete |
| Subsystem | retained scalar/selected-native C ABI and protocol boundary |
| Profiles | scalar, AVX-512, C, C++, Rust; CUDA and BTx/BTE excluded |
| Generation/provenance | Hand-maintained C++ boundary assembled from pinned upstream BLST entry points and CatBLST native implementations |
| Upstream references | Pinned BLST `54e6e55674722fc2797ebb4bbb71b26d881eb4b8`; symbol compatibility is tracked in `docs/audit/UPSTREAM_MANIFEST.md` |
| Related files | `bindings/blst_scalar.h`, `bindings/blst_types.h`, Rust bindings, selected/scalar point and field facades, `blst/sign/aggregate.cpp`, `blst/common/rb_tree.cpp`, scalar MSM, `test/C_api_tests.cpp`, `test/Aggregate_audit_tests.cpp`, `test/BindingHeader_tests.cpp` |

## Purpose and trust boundary

This file is the principal C ABI adapter for field-extension, group, MSM,
pairing, signature, key-generation, hash-to-curve, serialization, and
uniqueness operations. It accepts caller-owned objects, byte strings,
descriptor arrays, scratch storage, and opaque contexts. Wire bytes, messages,
DSTs, points awaiting validation, and all sizes and pointers are untrusted.
Scalars can be secret for constant-time multiplication, signing, and key
operations; MSM, FFT-like public scalar work, and explicit variable-time APIs
require public scalars.

The boundary must preserve the public layouts and symbols, perform no hidden
allocation in legacy `void` paths, prevent C++ exceptions from crossing C,
start C++ object lifetimes before typed access to caller storage, publish
fallible outputs according to the documented contract, and avoid conversions
on scalar-default hot paths. Selected AVX implementations are reached through
explicit boundary conversions only where their representation differs.

## Contracts and invariants

- Public C object representations are the exact `blst_*` types from
  `bindings/blst_types.h`; selected native C++ aliases and boundary static
  assertions pin compatible scalar layouts.
- Ordinary G1/G2/GT operations consume prime-order subgroup elements unless a
  declaration explicitly opts into the wider on-curve contract. The two
  inherited `blst_p1_unchecked_mult`/`blst_p2_unchecked_mult` functions are the
  explicit exact-scalar, on-curve exceptions. Release builds do not add hidden
  subgroup checks.
- `mult` is constant-time for canonical scalar inputs at fixed public `nbits`;
  supported noncanonical inputs are reduced modulo r with the documented
  canonical-versus-reduction timing distinction. `mult_vartime`, Pippenger,
  and weighted public operations require public scalars.
- Scalar bit lengths are at most 256. Pippenger window and fixed-table widths,
  scratch capacities, descriptor shapes, shifts, products, and pointer/length
  pairs are checked before dependent access.
- Compact pointer descriptors use a nonnull head followed by either explicit
  pointers or a null marker selecting contiguous continuation. Exact descriptor
  arrays require no trailing sentinel.
- Pippenger scratch is caller-owned `blst_limb_t` storage with the alignment
  and capacity returned by the corresponding query. At the C boundary every
  bucket lifetime is placement-started at a checked byte offset before typed
  access. Native already-live C++ bucket arrays retain their direct path.
- Empty ordinary G1/G2 Pippenger execution returns identity and its scratch
  query returns zero. Positive-count sizing and execution retain the upstream
  heuristic.
- `blst_miller_loop_n` is allocation-free. Scalar-default builds call the
  pointer-array implementation directly; selected AVX builds use bounded
  16-pair stack batches and publish output only after all batches complete.
- Opaque pairing and uniqueness contexts require storage sized/aligned by their
  query APIs. Initialization placement-starts the trivial fixed header and,
  for uniqueness, every variable-tail node. Reinitialization lawfully replaces
  these trivial objects without changing queried storage size.
- Serialization and uncompression establish on-curve representations, not
  implicit subgroup membership. Complete membership predicates test both the
  curve equation and subgroup relation. Advanced arithmetic remains narrow and
  requires its documented subgroup precondition.
- Exact output aliasing is supported only where declarations say so; partial
  overlap is not supported. Fallible functions leave output unchanged unless
  their declaration explicitly documents another result.
- No C++ exception may cross this boundary. The retained inventory found no
  throwing allocation or throwing helper on exported paths. Fallible work uses
  status, `bool`, null, bounded storage, or a documented serial/fail-closed
  fallback; legacy `void` paths are allocation-free/nonthrowing.
- Authoritative C and C++-callable declarations are in
  `bindings/blst_scalar.h`. The selected native C++ API is documented in the
  corresponding `blst/` facade headers. Raw Rust declarations are generated in
  `bindings/rust/src/bindings.rs`; safe wrappers and pairing ownership rules
  are documented in `bindings/rust/src/lib.rs`, `primitives.rs`, and
  `pairing.rs`. The declaration ledger in
  `docs/audit/PUBLIC_API_DOCUMENTATION.csv` records each family; no supported
  declaration surface implemented here is undocumented at the companion
  commit.

## Dependency and call-path map

Public C/C++ callers and safe Rust wrappers enter this translation unit. It
dispatches at compile time to selected `blst::G1Point`, `G2Point`, GT, map, and
field operations, or deliberately to scalar helpers when the public storage is
scalar-defined. Scalar MSM supplies affine conversion, fixed tables, and
Pippenger. Aggregate and key modules supply pairing contexts, aggregate verify,
and derivation. The uniqueness forwards terminate in `rb_tree.cpp`.

There is no runtime backend dispatcher or global thread pool. Scalar-default
paths use identical public/native types. AVX-default paths convert at the
declared boundary and retain native arithmetic internally; pure `cc` portable
emulation is not presented as the native public type.

## Upstream/specification comparison

The retained symbols and compact pointer conventions are compared against the
pinned upstream manifest. Intentional CatBLST divergences are: C/C++/Rust type
identity is expressed rather than assumed; subgroup and timing preconditions
are explicit; selected AVX paths remain compile-time; Miller batching is
allocation-free; invalid sizes fail closed; scratch/context object lifetimes
obey C++20; and CatBLST-only AVX target-group APIs use GT terminology. RFC 9380,
the BLS signature ciphersuite specifications, EIP-2333, and the public header
contracts govern hash/map, signing, verification, and derivation behavior.

## Checklist evidence

### 1. Identity and purpose

The complete file, all exported definitions, public declarations, selected
backend aliases, direct callees, tests, and history were inventoried. The C
boundary contains 253 retained definitions together with `exports.cpp`; every
definition is declared by the supported public header.

### 2. Contracts and invariants

Pointer/length pairs, scalar domains, subgroup/on-curve domains, timing,
aliasing, scratch shape, descriptor shape, opaque storage, error publication,
and output sizes were checked against C, C++-callable, and Rust surfaces. The
companion public-contract packet closes the missing declaration details.

### 3. Upstream and specification lineage

The upstream symbol inventory and behavior are pinned. CatBLST changes are
limited to documented contract, safety, selected-backend, naming, and
allocation-policy decisions; the arithmetic is owned and reviewed in its
lower-level files.

### 4. Correctness

Focused tests cover scalar/group boundaries, modulo-r multiplication,
serialization and membership separation, fixed-table and Pippenger layouts,
compact exact descriptor arrays, Miller products, aggregate flow, pairing
contexts, key derivation, hash/map vectors, and uniqueness behavior. Empty MSM
now agrees with execution by returning identity with zero required scratch.

### 5. Memory and lifetime safety

`AUD-186` repaired raw-limb Pippenger storage with a pointer-only workspace and
per-bucket lifetime start. Zeroing occurs while storage is raw, before objects
are started. `AUD-197` placement-starts pairing and uniqueness headers and
variable-tail nodes; the fake one-element trailing array is gone. Checked byte
offsets replace pointer arithmetic outside an actual array object. Temporary
AVX Miller batches are fixed-size stack arrays.

### 6. Undefined and implementation-defined behavior

Static assertions pin object size, alignment, trivial construction/destruction,
and selected type relationships used by the boundary. Shift counts, scalar
widths, scratch arithmetic, and descriptor counts are bounded. No C-style
whole-object type punning or indexing beyond separately created objects remains
in reviewed C scratch/context paths.

### 7. Side channels and cryptographic hygiene

Secret-capable multiplication and conditional operations retain fixed-schedule
arithmetic for canonical inputs. Public-variable-time algorithms are named or
documented as such. Release membership validation is not inserted into hot
paths. Selected AVX conversion occurs only at representation boundaries; no
runtime dispatch or hidden abstraction was added.

### 8. Randomness and secret lifecycle

Key generation and derivation use deterministic ciphersuite expansion and
erase their secret intermediates in their owning modules. Scalar preparation
wipes secret-derived temporaries. This file introduces no RNG and retains no
secret in global or heap storage.

### 9. API and ABI

Public symbols and signatures remain compatible with pinned BLST except for
previously approved CatBLST-only hard renames. `blst_pairing` and `blst_uniq`
have distinct opaque C tags without representation or symbol changes.
C++-only optional defaults are restored without affecting C declarations.
Rust raw results use ABI-safe integer newtypes and safe wrappers inherit the
documented C contracts.

### 10. Concurrency

All temporaries and scratch are caller-owned or per call. Pairing contexts are
mutable caller-owned objects and require external synchronization when shared.
There is no mutable process-global state, lazy initialization, or global
thread pool.

### 11. Resource exhaustion and denial of service

Size arithmetic is checked before multiplication, shifts, or access. Miller
uses bounded 16-pair batches. Pippenger and fixed-table APIs require queried
caller scratch and reject insufficient or invalid shapes. Lower-level batch
addition uses its reviewed bounded allocation/serial-fallback policy; no C ABI
exception or silent partial output results from allocation failure.

### 12. Portability and build behavior

The boundary uses fixed-width public types, requires 64-bit limbs, and is
compiled for scalar and selected AVX profiles. Strict C11/C++20 declaration
probes and MSVC/GCC/Clang syntax evidence cover the supported host surface.
Apple ARM64 remains a supported CI target. CUDA and BTx/BTE are intentionally
outside this review.

### 13. Performance-sensitive structure

Scalar-default operations remain direct calls. AVX-default operations stay
native after explicit boundary conversion. The Pippenger lifetime abstraction
inlines away: optimized MSVC retains the 0x108-byte core frame, reduces each
core instantiation by 16 bytes, and leaves G1/G2 tile wrappers at 0xE1/0x13D.
Pairing/uniqueness placement starts emit no helper calls or extra clearing.
No benchmark was warranted because inspected code generation did not change.

### 14. Tests and evidence

The exact lifetime packet's focused scalar Release MSM tests passed 12/12;
pairing/aggregate tests passed 9/9; uniqueness tests passed 3/3. Header and
public-contract tests passed 17/17 at the companion commit. Existing focused
boundary suites cover the remaining families. This documentation-only review
ran no build or test.

### Specialized additions

The C ABI exception inventory was performed across retained exports. CUDA,
BTx/BTE, benchmarks, and unsupported bindings were excluded by scope. Constant-
time and optimized assembly claims are recorded only where direct code-
generation evidence exists; algebraic correctness remains attributed to the
reviewed owning implementations.

## Findings

`AUD-042`, `AUD-048`, `AUD-069`, `AUD-074`, `AUD-080` through `AUD-082`,
`AUD-088`, `AUD-089`, `AUD-095`, `AUD-098`, `AUD-099`, `AUD-131`, `AUD-133`,
`AUD-134`, `AUD-137`, `AUD-139`, `AUD-149`, `AUD-152` through `AUD-154`,
`AUD-156`, `AUD-157`, `AUD-159`, `AUD-173`, and `AUD-185` through `AUD-194`,
plus `AUD-197`. In this final Wave 1 packet, `AUD-186`, `AUD-194`, and
`AUD-197` are validated.

## Validation

- Commit `02f5e5ff078b251ca3b5f86f4f83ea64af664bc7`: focused scalar Release
  MSM/Pippenger packet passed 12/12, including empty scratch and raw caller
  storage paths.
- Same commit: focused pairing/aggregate packet passed 9/9 and uniqueness
  packet passed 3/3.
- MSVC optimized-object inspection: no bucket construction loop, allocation,
  or out-of-line workspace call; core and tile frame/instruction shapes stated
  above were preserved.
- MSVC optimized-object inspection: pairing/uniqueness placement starts add no
  helper call or initialization beyond prior intended stores; ABI size queries
  are unchanged.
- Commit `aa4c7ed187f5a63496f0e94f9284108cbd46c7f1`: focused binding-header and
  contract packet passed 17/17, including distinct opaque tags and restored
  C++ defaults.
- No command was executed for this documentation-only record.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing retained-C-ABI baseline | Earlier issues mapped in `FINDINGS.md`; no unclosed Wave 1 blocker | Object lifetime, subgroup algebra, and constant-time code generation required hand review |
| Clang Static Analyzer / clang-tidy | Existing supported-host baselines | No additional blocker mapped to the final file | Cannot establish C++ lifetime of raw caller storage or cryptographic timing by itself |
| Clippy | Existing supported Rust binding runs | No safe-wrapper blocker after contract packet | Does not analyze the C++ implementation or raw FFI semantics |
| ASan/UBSan | Existing scalar baselines; not rerun for this packet | No mapped runtime report | Placement lifetime and selected AVX hardware evidence remain chiefly structural/focused |
| Fuzzing | Planned API and serialization campaign | Coverage gap, not a Wave 1 blocker | Full malformed descriptor/wire corpus remains milestone work |

## Residual risk and decision

The file is complete at the recorded implementation commit with public
declarations reviewed at the companion commit. Residual release evidence is
the planned larger sanitizer/fuzz/coverage matrix, native AVX-512 cross-host
execution, and continuing Apple CI. Experimental CUDA and frozen BTx/BTE are
not covered. These gaps do not reopen the exact scalar C boundary because its
contracts, exception paths, object lifetimes, resource policy, and focused
code generation have been reviewed directly.
