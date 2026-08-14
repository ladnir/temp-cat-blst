# File Review: `blst/scalar/fields/Fp12.cpp`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/scalar/fields/Fp12.cpp` |
| Roadmap phase / queue | Phase 14 scalar extension fields and GT boundary |
| Reviewed commit | `2333650dcd3017f34e169c79bc8dc86e27d7305e` |
| Reviewer | Codex audit team |
| Date | 2026-08-13 |
| Disposition | in-scope production |
| Review state | complete |
| Subsystem | scalar Fp6/Fp12 arithmetic, final exponentiation, GT GLV exponentiation |
| Profiles | scalar, C, C++, Rust; experimental CUDA consumers excluded from release evidence |
| Generation/provenance | Hand-maintained C++ refactor of upstream tower/pairing code plus local Hypercat GT/GLV work |
| Upstream references | `54e6e556:src/fp12_tower.c`, `src/pairing.c`; local GT/GLV lineage begins before scalar refactor `0eba90d` |
| Related files | `bindings/blst_types.h`, `bindings/blst_scalar.h`, `blst/scalar/curve/GT.h`, `blst/c_api.cpp`, Fp2/Fp, scalar Fp12/GT tests |

## Purpose and trust boundary

This file implements the Fp6/Fp12 tower, sparse Miller multiplication,
cyclotomic squaring, inverse, Frobenius, final exponentiation, serialization,
and two- and four-dimensional GT exponentiation. Arithmetic values can cross
the public C/C++/Rust boundary. GT APIs distinguish fixed-schedule secret-
capable multiplication from explicitly variable-time public-scalar operations.
All storage is fixed stack/static; the file has no heap allocation, ordinary
exception source, callback, or mutable global state.

## Contracts and invariants

- General Fp12 arithmetic accepts reduced Montgomery coefficients.
- Cyclotomic square/exp requires the cyclotomic subgroup; GLV additionally
  requires the order-r GT subgroup.
- General multiply, square, sparse multiply, inverse, Frobenius, and final-exp
  support exact whole-object aliases where exposed.
- Frobenius is periodic with degree 12; the inherited coefficient helper only
  directly represents degrees one through three.
- GLV schedules must be compiler-produced and within their fixed capacities.
- Raw 576-byte input currently reduces coefficients rather than validating a
  canonical encoding.

## Batch B widened-arithmetic delta

The widened Fp2 helpers require each component below `p * 2^384`. Debug builds
now assert this input bound. The `u+1` helper also asserts its existing
non-alias contract. Release builds erase these checks. A focused MSVC Debug
compile accepts the resulting Fp12 path, while the byte-level Release oracle
checks the underlying wide add/subtract operations and both exact aliases.

## Dependency and call-path map

The tower composes Fp2/Fp primitives and selected widened scalar assembly.
Pairing/final verification consumes the inherited arithmetic. Native GT and
the public `blst_gt_*` wrappers consume local GLV paths. The Fp12 C wrappers are
allocation-free and nonthrowing by construction; Rust reaches the same C ABI.
Scalar and AVX GT implementations are independent and require differential
evidence when semantics change.

## Upstream/specification comparison

Core tower formulas, constants, sparse multiplication, Frobenius coefficients,
cyclotomic square, group predicates, and final exponentiation retain pinned
upstream structure. Upstream explicitly restricted the exposed Frobenius degree
to one through three. Local C++ documentation omitted that restriction, and the
local GT/GLV exponentiation/schedule/serialization APIs have no upstream BLST
equivalent. Their subgroup and timing contracts therefore require explicit
CatBLST design decisions.

## Checklist evidence

### 1. Identity and purpose

All 1,990 lines, declarations, direct supported callers, tests, and local
history were inspected at the recorded commit.

### 2. Contracts and invariants

Zero/one, aliasing, exponent width, Frobenius degree, subgroup, schedule,
serialization, and inverse behavior were traced. Violations are recorded in
`AUD-079` through `AUD-086`.

### 3. Upstream and specification lineage

Core routines were compared with pinned `fp12_tower.c` and `pairing.c`. Local
GLV code was separated from inherited arithmetic. The file itself lacks the
required license/provenance header (`AUD-087`).

### 4. Correctness

Inherited fixed-width tower formulas and alias-ordering remain structurally
sound. General cyclotomic exponentiation incorrectly selects an order-r-only
GLV identity, and the public Frobenius degree is unchecked. The declared
`glv_split_4` operation has no definition.

### 5. Memory and lifetime safety

Typed fixed stack storage replaces the prior identity union. No heap or retained
pointer exists. Unchecked Frobenius and schedule indices are concrete
out-of-bounds reads. GLV4 uses roughly 9 KiB for its 16-entry Fp12 table before
additional temporaries; this is accepted for supported hosts but must not grow
without stack/code-generation evidence.

### 6. Undefined and implementation-defined behavior

Little-endian and 64-bit limbs are globally required. Prior wrapped signed-digit
conversion uses `bit_cast`. Remaining narrowing, signedness, and shift warnings
are listed in `AUD-087`; exact signed-component bounds need regression evidence.

### 7. Side channels and cryptographic hygiene

Tower arithmetic has fixed structure. Euclidean inversion remains
variable-time. GT exponentiation now has a fixed-schedule safe default and an
explicitly named variable-time GLV schedule for public scalars. Optimized
code-generation inspection and focused differential tests are recorded under
`AUD-082`/`AUD-088`.

### 8. Randomness and secret lifecycle

N/A: no entropy. The variable-time path requires public scalars. The
constant-time path wipes its scalar decomposition and magnitude temporaries
after use.

### 9. API and ABI

C wrappers use exact public representations and have no allocation/exception
path. `blst_fp12_frobenius_map` checks its deliberately narrow 1-3 degree
domain before calling the native primitive. Native C++ adds schedule and
deserialize APIs with incomplete safety contracts. Future
symbol-prefix configuration does not affect the internal hot-path design.

### 10. Concurrency

All constants are immutable and there is no lazy cache, thread pool, or shared
mutation. Calls are reentrant subject to caller-owned output alias contracts.

### 11. Resource exhaustion and denial of service

Supported C GT scalar length is documented at at most 256 bits but not enforced.
General C++ exponent work is caller-sized. Invalid schedule lengths can extend
work and memory reads beyond fixed storage until remediated.

### 12. Portability and build behavior

MSVC and Clang parse the selected host file. Strict MSVC reports four local
conversion/shift warnings; strict Clang reports the corresponding conversions
and signedness cases. Experimental CUDA is not release evidence.

### 13. Performance-sensitive structure

Widened Fp2x2/Fp6 arithmetic, sparse multiplication, cyclotomic square ordering,
fixed addition chains, and native calls are accepted hot structure. Domain and
bounds fixes must keep the one-to-three Frobenius path direct, validate schedules
once outside inner loops, allocate only on exact alias paths, and retain GT-only
GLV selection where approved.

### 14. Tests and evidence

Current tests cover inherited arithmetic, Frobenius boundary degrees, and both
reduced and full-width GT/GLV scalars. They do not cover arbitrary cyclotomic
non-GT bases, malformed schedules, noncanonical input, or non-256-bit native
in-place GT exponentiation.

### Specialized additions

Scalar/AVX GT differential tests are required for approved exponent changes.
Assembly generators are unchanged by this review. CUDA tests are explicitly not
a release gate.

## Findings

`AUD-079`, `AUD-080`, `AUD-081`, `AUD-082`, `AUD-083`, `AUD-085`, `AUD-086`,
`AUD-087`, `AUD-088`, `AUD-198`, `AUD-204`, `AUD-205`, `AUD-206`, and
`AUD-207`.

## Validation

- Full source, declaration, caller, test, upstream, and history inspection.
- `cl /std:c++latest /W4 /WX /Zs ... Fp12.cpp`: four local warnings, recorded
  in `AUD-087`; no object/build work performed.
- `clang++ -std=c++20 -Wall -Wextra -Wpedantic -Wconversion
  -Wsign-conversion -fsyntax-only ... Fp12.cpp`: parses; local conversions map
  to `AUD-087`, shared-header warnings remain in their own reviews.
- Structural searches prove missing `glv_split_4`, unused functions/macros,
  dead branches, and absent negative/boundary tests.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing C/C++ baseline | No separate new alert | Subgroup algebra, timing, and API preconditions require manual review |
| Clang Static Analyzer / clang-tidy | Not rerun | Pending remediation | Strict compiler diagnostics only in this packet |
| Clippy | N/A | N/A | Rust reaches the C implementation |
| ASan/UBSan | Prior scalar baseline | Does not exercise invalid degree/schedule | Focused negative cases required |
| Fuzzing | Not yet present | Coverage gap | Schedule, serialization, and differential Fp12 targets are candidates |

## Residual risk and decision

### Post-review remediation

The following decisions supersede the open-state descriptions above:

- `AUD-079` is fixed. Zero-bit exponentiation returns one, reverse iteration is
  `size_t`-safe, and exact Fp12 aliasing uses a non-inlined preservation helper.
  Optimized MSVC leaves the normal path with a 40-byte frame and no base copy.
- `AUD-080` is fixed with a checked C/Rust boundary that returns `false` and
  preserves output outside powers 1-3. The internal primitive keeps its direct
  fixed-table structure and explicit precondition.
- `AUD-081` documents the historically named `cyclotomic_exp` operation as
  GT-specific. The GLV implementation accepts the complete 256-bit unsigned
  scalar domain modulo the GT order; direct basic-exponent comparisons pass at
  r, r+1, 2r-1, 2r, and 2^256-1 in scalar and native AVX-512/IFMA execution.
  Both splitters assert the exact final quotient bound of two in Debug.
- `AUD-082`/`AUD-088` make unsuffixed GT multiplication fixed-schedule and
  retain an explicit optimized `mult_vartime`. CT uses four unsigned radix-z
  digits, conjugated odd Frobenius bases, a final 16-entry owning table, and a
  fixed 64-round full-scan executor. Scalar, portable-AVX, native AVX, and Rust
  differential checks pass. GCC 15.2 optimized inspection finds no gathers,
  allocation references, or out-of-line selected-vector wrappers; the native
  GT frame is 9,984 bytes. Sequential Peach measurements and the scalar/AVX
  variable-time splitter decisions are recorded in `PHASE_14.md`.
- `AUD-083` is accepted as a low-level precondition: schedules are mutable,
  variable-time trusted bytecode, and only unmodified compiler output is valid.
- `AUD-085` is fixed. `Fp12::from_bendian` rejects every coefficient at least
  p and publishes only a complete canonical Montgomery value. Native checked
  GT decoding adds membership validation; `deserialize_unchecked` deliberately
  defers only membership for aggregate validation.
- `AUD-086` is fixed by removing the obsolete undefined duplicate declaration.
- `AUD-087` is fixed: provenance, macro names/leakage, inactive branches,
  commented code, placeholder tests/comments, and local warnings were cleaned.

Focused MSVC validation passes zero-bit Fp2/Fp12 exponentiation, exact Fp12
aliasing, and canonical Fp12 serialization. Scalar Debug/Release and selected
portable-AVX Debug pass 5/5 GT decoder cases covering p-1, p in every
coefficient position, canonical non-GT input, failure atomicity, and deferred
aggregate checking. `cl /W4 /WX /Zs` passed the earlier implementation packet.
No broad suite or decoder benchmark ran.

The timing decision is resolved by constant-time default and named variable-
time GT multiplication (`AUD-082`/`AUD-088`). The checked versus deferred
deserialization contract (`AUD-085`) is also resolved. Complete source and
caller review found no remaining implementation decision. Supported run
`31665323681` passes the containing revision on Apple Silicon, Linux, and
Windows. The file review is complete; broader sanitizer, oracle, and fuzz
campaigns remain subsystem release gates.

## Wave 2 exact-contract follow-up

The implementation and all direct supported callers were reread at `314c186`.
No arithmetic change was required. The invalid-domain cyclotomic-square test
now final-exponentiates its seed and compares the specialized result with
general squaring. Predicate, exponent, schedule, and sparse-operand declaration
contracts now match the implementation. Broader sanitizer and independent
oracle evidence remains.
