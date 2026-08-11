# File Review: `blst/common/timer.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/common/timer.h` |
| Roadmap phase / queue | Phase 13/17, optional FFT and middle-product diagnostics |
| Reviewed commit | Packet 17C working tree based on `a0fe5c801090190463876bd4d2add18e22fa1235`; final commit pending |
| Reviewer | Codex audit team |
| Date | 2026-08-03 |
| Disposition | `in_scope_production` |
| Subsystem | Caller-owned optional timing diagnostics |
| Profiles | scalar and AVX C++; C APIs pass a null sink |
| Generation/provenance | Locally composed C++ header |
| Upstream references | `not_applicable: local FFT/middle-product instrumentation` |
| Related files | `blst/fft/ctFft.h`, `blst/fft/middle_product.h`, `test/Timer_tests.cpp`, timer audit probes |

## Purpose and trust boundary

This header optionally accumulates named durations for native C++ FFT and
middle-product diagnostics. Product calls commonly pass `nullptr`; that path
must be almost free. Diagnostics must never throw, terminate crypto work, or
introduce shared global state.

## Contracts and invariants

- A null `TimingStats*` performs no clock read or bookkeeping.
- Construction, destruction, simple/composed `add`, `get_ns`, and `clear` are nonthrowing.
- Allocation failure drops the new diagnostic sample without changing product
  output or propagating an exception.
- Keys referenced by `ScopedTimer` must outlive that scope.
- `TimingStats` is caller-owned, non-copyable, non-movable, and not internally
  synchronized.
- Enabled timing is diagnostic and not a constant-time or benchmark facility.

## Dependency and call-path map

FFT and middle-product functions create stack `ScopedTimer` objects and accept
an optional caller-owned sink. Retained C APIs pass null. No global pool,
registry, background task, or public C timing state exists.

## Upstream/specification comparison

N/A: this is local instrumentation. The review compares behavior to the local
security/performance policy rather than an upstream implementation.

## Checklist evidence

### 1. Identity and purpose

Every construction site and sink flow in FFT/middle-product code was searched.

### 2. Contracts and invariants

Null, exception, ownership, key lifetime, and synchronization rules are listed
above and encoded in the types.

### 3. Upstream and specification lineage

N/A with rationale above.

### 4. Correctness

Repeated keys accumulate nanoseconds, absent keys return zero, and `clear`
returns to an empty lazy state. Transparent lookup avoids diagnostic read
allocation.

### 5. Memory and lifetime safety

The optional map is constructed lazily and destroyed by its owner. The timer
stores only a non-owning key view for its lexical lifetime. Copy/move is
disabled to prevent duplicate final recording.

### 6. Undefined and implementation-defined behavior

Steady-clock arithmetic uses standard chrono types. Signed duration is cast to
the existing unsigned statistic after a forward steady-clock interval.

### 7. Side channels and cryptographic hygiene

Enabled timing intentionally observes execution duration and must not be
exposed across an adversarial secret boundary. Normal C/product calls pass
null; the disabled optimized path is empty.

### 8. Randomness and secret lifecycle

N/A: no entropy or secret storage. Key strings describe public code regions.

### 9. API and ABI

Native C++ diagnostic types are internal and have no stable public ABI. The C
API receives no timer object.

### 10. Concurrency

There is no global state or lock. A sink may be used only according to ordinary
caller-owned synchronization; concurrent mutation of one sink is unsupported.

### 11. Resource exhaustion and denial of service

Enabled unique keys can allocate, but failure is contained and callers choose
the instrumentation. Production null calls allocate nothing.

### 12. Portability and build behavior

Standard C++ chrono/container code builds under MSVC 19.50. No CUDA consumer is
part of the supported timing surface.

### 13. Performance-sensitive structure

MSVC `/O2 /Ob3` compiles a lexical null timer to exactly one `ret`, with no
clock call, stack traffic, spill, allocation, or exception path. Enabled timing
retains one start/end clock pair and map update.

### 14. Tests and evidence

No-throw traits, accumulation, missing key, clear, null scope, optimized code
generation, and forced allocation failure are covered by dedicated probes.

### Specialized additions

FFT mathematics is unchanged; this review covers instrumentation only.

## Findings

`AUD-065` and `AUD-296` - validated.

## Validation

- `Timer.*`: 1/1 passed, including composed prefix/suffix keys.
- MSVC 19.50 `/O2 /Ob3 /W4 /WX` null probe: exact `ret`.
- Forced global-allocation failure probe: exited zero, preserved previous
  simple and composed samples, and propagated no exception.
- Single-worker MSVC Debug target build: passed.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Phase 12 database | `AUD-065` hand audit | Zero-cost behavior requires optimized object inspection |
| Clang Static Analyzer / clang-tidy | Retained profiles | No additional mapped diagnostic | Does not prove allocation-failure containment at destruction |
| Clippy | N/A | N/A | Native C++ diagnostics only |
| ASan/UBSan | Prior baselines | No mapped report | Allocation injection and codegen provide direct evidence |
| Fuzzing | N/A | N/A | No untrusted grammar; finite state covered directly |

## Residual risk and decision

Complete at the recorded commit. Enabled timing can perturb and reveal timing,
is not synchronized, and can wrap a 64-bit accumulated counter; those are
documented diagnostic limitations. Release callers should keep the sink null
unless explicitly collecting local evidence.
