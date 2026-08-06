# File Review: `blst/scalar/curve/point.cpp`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/scalar/curve/point.cpp` |
| Roadmap phase / queue | Phase 15 scalar curves and protocols, first packet |
| Reviewed commit | `07d79a48cc7d82560f5194f76a137e4ee2527d06`; localized fixed-window follow-up `b9855299349fa083c455d45d8cbbb30f2a617e03` |
| Reviewer | Codex audit team |
| Date | 2026-08-04 |
| Disposition | in-scope production |
| Subsystem | scalar G1/G2 constants, arithmetic, serialization, subgroup relations, signing helpers |
| Profiles | scalar, C, C++, generated/safe Rust callers; experimental CUDA excluded from release gate |
| Generation/provenance | Hand-maintained C++ translation of upstream BLST point formulas plus local GLV/GLS work |
| Upstream references | pinned BLST `src/ec_ops.h`, `src/e1.c`, `src/e2.c`; CatBLST pre-refactor macro history |
| Related files | `point.h`, scalar fields, `GT.cpp`, `c_api.cpp`, map/sign/MSM callers, focused point tests |

## Purpose and trust boundary

This file defines scalar point constants, complete and incomplete group
formulas, normalization/encoding/decoding, membership relations,
endomorphisms, scalar multiplication, signing helpers, and XYZZ conversions.
It receives untrusted bytes only through decoders; arithmetic inputs otherwise
rely on the documented caller-established curve/subgroup contracts.

The file has no heap allocation, exception source, callback, mutable global
state, or lazy initialization. Generator constants are immutable and typed.

## Contracts and invariants

- Field coordinates are reduced Montgomery Fp/Fp2 values.
- Infinity and subgroup/timing/alias rules match the reviewed header record.
- Incomplete formulas are selected only through narrow `add`; `dadd` handles
  equal/opposite inputs.
- Raw decoders ignore the top three x bits and return encoded errors where
  inherited; high-level native decoders publish only on success.
- Generator functions return references to live objects of the declared type.
- Projective serialization uses typed x/y references after optional Jacobian
  normalization; it does not form an affine reference to projective storage.

## Dependency and call-path map

Fp/Fp2 operations and scalar assembly implement field arithmetic. C wrappers
convert exact boundary objects and copy local results. Scalar GT/pairing uses
typed affine generators. Hash-to-curve, signing, aggregation, and MSM call the
point primitives. `bulk_addition.h` was traced as a caller and yielded the
separate open `AUD-109` scratch-lifetime finding.

## Upstream/specification comparison

EFD complete/incomplete formulas and generator constants retain inherited
structure. The old C macro supplied `get_BLST_Rx()->p2` for the G2 identity;
the C++ translation collapsed that to `Fp::one()`, causing `AUD-104`. Prefix
casts were inherited layout techniques but are not valid live-object creation
in standard C++ (`AUD-105`).

## Checklist evidence

### 1. Identity and purpose

All 2,061 source lines, declarations, supported callers, tests, and relevant
history were inspected.

### 2. Contracts and invariants

Arithmetic domains, encoding flags, output publication, exact aliasing,
constant-/variable-time boundaries, and generator lifetimes were traced.

### 3. Upstream and specification lineage

The pre-refactor macro/union operand selection and later typed refactors were
compared to explain the G2 regression rather than merely patch its symptom.

### 4. Correctness

`Field::one()` now supplies the correct G1 or G2 z identity. Typed generator
coordinates match their Jacobian generators. Serialization round-trips and
complete equal/opposite behavior pass focused tests.

### 5. Memory and lifetime safety

G2 mixed addition no longer reads past Fp storage. Generator/serialization and
Miller fallback paths no longer create nonexistent affine objects. Native G1
decoder failure does not partially publish output.

### 6. Undefined and implementation-defined behavior

The reviewed file has no remaining Point-to-Affine reference cast. Raw byte
operations use existing audited vector/byte primitives. Cross-file scratch
prefix casting is `AUD-109`.

### 7. Side channels and cryptographic hygiene

Add/dadd/doubling and conditional operations retain branchless structure.
Conversion/serialization/membership branches operate on public points.
Canonical multiplication timing and scalar erasure are covered by
`AUD-088`/`AUD-090`.

### 8. Randomness and secret lifecycle

N/A for randomness. Signing consumes validated secret-key storage through the
reviewed fixed-schedule multiplication path; local scalar erasure is recorded
in the multiplication findings.

### 9. API and ABI

No C signature changed. The undefined internal `mult_w5` declaration was
removed. C/Rust raw bindings retain inherited names and exact layouts.

### 10. Concurrency

All globals are immutable constants. Calls are reentrant with caller-owned
outputs and no shared cache or thread pool.

### 11. Resource exhaustion and denial of service

All operations here use fixed-size storage and bounded scalar widths. Encoding
input sizes are fixed by the declaration contracts.

### 12. Portability and build behavior

MSVC Release/Debug and GCC ASan compile and execute the focused target. Typed
CUDA affine constants were kept structurally consistent, but CUDA remains
experimental and was not built.

### 13. Performance-sensitive structure

The arithmetic fix changes only a constant address. Serialization retains one
normalization Point and passes x/y references directly; optimized MSVC shows no
affine copy or extra buffer. No benchmark was run because code shape did not
materially change.

### 14. Tests and evidence

Focused tests cover G1/G2 typed generators, serialization/compression parity,
G2 mixed identity and exact aliasing, complete equal/opposite addition, native
decoder failure preservation, and valid publication.

### Specialized additions

No assembly generator changed. Apple scalar evidence remains inherited from
the supported Apple CI gate; this packet did not trigger a new CI run.

## Findings

`AUD-024`, `AUD-035`, `AUD-070`, `AUD-088`, `AUD-089`, `AUD-093`, `AUD-096`,
`AUD-104`, `AUD-105`, `AUD-107`, `AUD-108`, `AUD-127`, `AUD-130`; caller
finding `AUD-109`.

## Validation

- MSVC 19.50 Release scalar focused build/execution: 6/6 passed.
- MSVC 19.50 Debug scalar focused build/execution: 6/6 passed.
- Linux GCC 13 ASan/LSan Debug focused build/execution: 6/6 passed.
- Optimized MSVC disassembly: one Point normalization frame; coordinate helper
  receives direct x/y addresses; no affine copy/second point buffer.
- Structural and documentation/group-contract checks: passed.
- The localized `ec_mult.h` follow-up consumes the precomputed status and
  removes wipes only from the public variable-time wNAF path; focused Release
  execution and optimized inspection are recorded with that packet.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing C/C++ baseline | No separate point alert | Object lifetime, formulas, and timing required manual trace |
| Clang Static Analyzer / clang-tidy | Existing baseline | No new blocker | CUDA excluded; algebraic preconditions not modeled |
| Clippy | Existing supported profiles | No direct implementation | Raw Rust reaches C ABI; safe containing file pending |
| ASan/UBSan | GCC 13 ASan/LSan Debug | 6/6 passed, including AUD-104 | UBSan not rerun |
| Fuzzing | Not yet present | Coverage gap | Compressed/uncompressed decoding is a priority target |

## Residual risk and decision

The original full review and localized fixed-window follow-up are complete at
their recorded commits. The later `AUD-094`/`AUD-160` declaration pass defines
successful point decode as on-curve only and reserves optimized/advanced
operations for subgroup inputs. Focused non-subgroup and decoder-flow evidence
remains; no arithmetic changed here. GT decoding policy (`AUD-085`) was
later resolved by the GT/Fp12 review: the default decoder validates membership,
while the explicitly unchecked decoder defers only that check.
