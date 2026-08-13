# File Review: `blst/avx/curve/ec.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/avx/curve/ec.h` |
| Roadmap phase / queue | Phase 16 AVX curve and MSM closure |
| Reviewed commit | `5b08edc8f890ef1cbf4c6443162ec4da985879ee` |
| Reviewer | Codex audit team |
| Date | 2026-08-13 |
| Disposition | `in_scope_production` |
| Subsystem | AVX-512 curve arithmetic |
| Profiles | native AVX-512 and portable developer Debug |
| Generation/provenance | VROOM-derived MIT point carriers and doubling; locally added a4 schedule |
| Upstream references | SimonLangowski/VROOM `d6ecd534c54a41702ea80ae75b52250b6d51b372` |
| Related files | `matrixnok_curve.h`, `matrixnok_g1.h`, `matrixnok_g2.h`, `point.h`, `scalar_mult.h`, AVX corner-case tests |

## Purpose and trust boundary

The file owns fixed-storage homogeneous and affine point carriers plus the
low-level schedules still selected by the MatrixNoK G1/G2 facades. It performs
no point validation. Callers supply points on the stated curve and the ring
whose reduction and bound contracts match their coordinates.

Projective infinity is semantic `Z = 0`; selected construction emits
`(0,1,0)`. Affine infinity is the exact `(0,0)` sentinel. Raw RNS zero helpers
therefore remain explicitly representation-level observations.

## Contracts and invariants

- `PointDouble` accepts any on-curve projective point, including infinity.
  Subgroup membership is unnecessary. The ring must implement the curve's
  `3b` multiplication and accepted MatrixNoK bounds.
- `PointDaddA` accepts two points on the same curve. Its `a4` argument is that
  curve's coefficient, not a tuning parameter. The owning facade performs the
  identity selection.
- `ProjectiveEqual` requires the supplied semantic-zero predicate and compares
  homogeneous cross-products after intentional RNS normalization.
- The retained schedules allocate nothing, throw nothing, and contain no
  mutable global state or runtime dispatch.
- Point-data choices use mask selection. Choosing the public `a4` schedule is
  allowed to branch before entering this file.
- Exact output/input aliasing is safe through value-returning schedules.

## Dependency and call-path map

`matrixnok_g1.h` and `matrixnok_g2.h` select `PointDouble`. `point.h` selects
`PointDaddA`, exact identity helpers, and predicate-based equality. Pairing and
scalar-multiplication code use the carriers but select their own arithmetic
facades. The deleted generic G1/G2 engine had no retained caller.

## Upstream/specification comparison

The carriers and `PointDouble` descend from VROOM commit `d6ecd534` under MIT.
`PointDouble` implements Algorithm 9 of ePrint 2015/1060 and retains its exact
batched schedule. Category Labs added `PointDaddA` during AVX integration. Its
correctness is covered by scalar differential tests, not by the complete-add
Lean theorem. `THIRD_PARTY_NOTICES` now reproduces the required VROOM notice.

## Checklist evidence

### 1. Identity and purpose

Every retained type and operation has a live selected caller. Zero-caller
generic formulas, facades, defaults, and replacement helpers were removed.

### 2. Contracts and invariants

Coordinate representations, point domains, `a4`, allocation, exception,
timing, and alias contracts are stated at the owners.

### 3. Upstream and specification lineage

The header pins VROOM `d6ecd534`, distinguishes local `PointDaddA`, uses the
MIT SPDX identifier, and retains Category Labs attribution.

### 4. Correctness

The accepted Lean row continues to bind selected complete addition and public
doubling. New scalar differentials cover both fields and all important a4
addition shapes.

### 5. Memory and lifetime safety

The file contains owning value types only. It performs no pointer arithmetic,
allocation, placement construction, or type punning.

### 6. Undefined and implementation-defined behavior

No integer shift, unbounded index, representation cast, or signed-overflow
operation is present. Ring arithmetic owns its separate machine bounds.

### 7. Side channels and cryptographic hygiene

Schedules are fixed. Point selection delegates to selected mask operations.
Compiler lowering and microarchitectural leakage remain P0 evidence boundaries.

### 8. Randomness and secret lifecycle

N/A: the file consumes no randomness and owns no standalone scalar buffer.

### 9. API and ABI

These are internal native C++ templates. Public types and C symbols are owned
by `point.h` and `bindings/blst_avx.h`.

### 10. Concurrency

All state is caller-owned. No globals, locks, caches, or shared mutation exist.

### 11. Resource exhaustion and denial of service

All storage and work are fixed by the selected formulas.

### 12. Portability and build behavior

Portable-Debug MSVC compiles every affected owner. Native execution remains a
separate AVX-512/IFMA release profile.

### 13. Performance-sensitive structure

Selected arithmetic bodies are unchanged. Deletion adds no abstraction,
conversion, copy, branch, or call to a hot path.

### 14. Tests and evidence

A focused executable passed all eight AVX add-corner cases. Structural searches
found no removed curve facade or implicit scalar-multiplication default.

### Specialized additions

`AVX-P4-003` is rebound at the reviewed hash. The mapped public-doubling body
is unchanged; `PointDaddA` remains explicitly differential-only.

## Findings

`AUD-170`, `AUD-225`, `AUD-228`, `AUD-241`, `AUD-387`, `AUD-388`.

## Validation

- Focused portable-Debug compilation of the AVX C API and three test objects.
- Focused execution: 37/37 curve, MSM, and C-boundary tests passed; 8/8 belong
  to `AvxAddCornerCases`.
- `validate_avx_correspondence.py`: 21 rows and 61 hashes passed.
- Structural searches and `git diff --check` passed.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing C/C++ baseline | no mapped defect | Template arithmetic and constant time need direct evidence |
| Clang analyzer / clang-tidy | Existing AVX triage | no mapped defect | Does not prove ring formulas |
| Clippy | N/A | N/A | C++ owner |
| ASan/UBSan | Not rerun | deferred | No lifetime or indexing code changed |
| Fuzzing | Not run | deferred | Differential tests cover selected shapes, not a corpus |

## Residual risk and decision

Native compiler and ISA behavior, microarchitectural constant time, and a
machine proof of `PointDaddA` remain outside this review. The complete source
review is closed because every retained operation, caller, and contract has
been traced and focused differentials pass.
