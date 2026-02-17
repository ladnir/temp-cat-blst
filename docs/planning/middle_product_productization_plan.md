# Middle Product Productization Plan

## Goal
Productize the middle-product workflow currently exercised in benchmarks into reusable library code under `blst/fft/`, with two callable algorithm variants:
1. Basic FFT path (full convolution path)
2. Specialized middle-product path (reduced work path)

The runtime winner depends on input size, so both implementations must be exposed and test-covered.

## Scope
In scope:
- Add public internal API in `blst/fft/` for both middle-product algorithms.
- Refactor benchmark logic into library-quality functions (no benchmark-only coupling).
- Add correctness tests with a deterministic `O(n^2)` reference for the middle window.
- Document CPU-only status and CUDA follow-up plan.

Out of scope for first milestone:
- Selecting a default algorithm automatically at runtime.
- Production CUDA implementation parity.
- Large API redesign outside FFT/middle-product surface.

## Target API (Phase 1 design)
Add a new header/source pair in `blst/fft/` (exact filenames to finalize during implementation), with an API shape similar to:

```cpp
enum class MiddleProductAlgorithm {
    BasicFFT,
    SpecializedFFT,
};

struct MiddleProductOptions {
    MiddleProductAlgorithm algorithm;
    bool b_is_pretransformed;  // offline FFT(B) already available
};

// Computes middle n coefficients for A in G1 and B in G2.
// Returns output in GT.
void middle_product(
    std::span<const G1Point> a,
    std::span<const G2Point> b_or_bhat,
    size_t n,
    std::span<GT> out,
    const MiddleProductOptions& opts);
```

Notes:
- `b_or_bhat` supports either raw B or pretransformed FFT(B), controlled by options.
- Keep helper APIs internal if only used by tests/benchmarks.
- Preserve existing benchmark compatibility by adding thin wrappers where needed.

## Implementation Phases

### Phase 0: Baseline and Extraction Map
- Identify benchmark functions to extract from `benchmark/FFT_bench.cpp`:
  - `bench_middle_product_basic`
  - `bench_middle_product_fft`
  - Supporting helpers used by both
- Define ownership boundaries between benchmark harness and reusable core logic.
- Record required precomputed artifacts (roots, inverse roots, primitive root, scratch buffers).

### Phase 1: Core CPU Library Implementation
- Add middle-product core implementation under `blst/fft/`.
- Implement both variants behind one shared interface:
  - Basic FFT full pipeline
  - Specialized middle-product pipeline
- Explicitly model offline/online flow for B precompute:
  - Provide helper to precompute FFT(B) on CPU.
  - Support online call consuming precomputed FFT(B).
- Ensure no benchmark dependency remains in core files.

### Phase 2: Benchmark Rewire
- Update benchmark code to call new core APIs only.
- Preserve benchmark names for continuity of historical comparisons.
- Keep benchmark parameter sets intentionally small for dev speed, with filter-driven full sweeps later.

### Phase 3: Correctness Tests
Add tests under `test/` focused on correctness of the middle window:

1. Reference comparison (basic path)
- Build scalar coefficient inputs over `Fr`.
- Compute expected middle coefficients using `O(n^2)` field multiplication.
- Project coefficients to group setting and compare resulting GT outputs.

2. Reference comparison (specialized path)
- Same reference oracle as above.
- Validate specialized algorithm output matches oracle exactly.

3. Cross-algorithm equivalence
- For same random seeds and sizes, verify BasicFFT == SpecializedFFT.

4. Pretransformed-B path
- Validate results are identical when using raw B vs precomputed FFT(B).

5. Edge and boundary cases
- Small n (1, 2, small powers of two, non-powers)
- n values that change FFT size bucket
- Deterministic seed-based randomized regression set

Test guidance:
- Use deterministic vectors and fixed seeds.
- Keep tests CPU-only for this milestone.
- Avoid benchmark timing assertions in correctness tests.

### Phase 4: Docs and Adoption
- Add short usage notes in FFT docs or README section:
  - How to choose algorithm explicitly
  - How to precompute FFT(B) offline
  - CPU-only support status
- Add TODO markers and tracking item for CUDA enablement.

## CUDA Status and Forward Plan
Current status:
- Productized middle-product API target is CPU-first.
- CUDA benchmark and support paths are lagging and not part of milestone acceptance.

CUDA follow-up backlog:
- Add `#ifdef BLST_CUDA` extension hooks in new API surface.
- Mirror BasicFFT and SpecializedFFT kernels for GPU path.
- Reuse same correctness tests with CPU-vs-GPU parity checks once stable.
- Add benchmark labels clearly indicating CPU vs CUDA implementation.

## Acceptance Criteria
- Both algorithm variants callable from library code (not only benchmark code).
- Existing benchmark targets compile and run against new API.
- New tests pass and validate middle window correctness against `O(n^2)` oracle.
- Offline FFT(B) path is represented in API and tested.
- CPU-only limitation and CUDA roadmap documented.

## Risks and Mitigations
- Risk: Scratch-buffer ownership/performance regressions.
  - Mitigation: Explicit scratch object or controlled internal allocation strategy.
- Risk: Ambiguity between raw B and FFT(B) call paths.
  - Mitigation: Strong option typing and separate helper for pretransform.
- Risk: Algorithm mismatch on edge sizes.
  - Mitigation: Boundary-focused test matrix and cross-algorithm equivalence tests.

## Proposed Execution Order
1. Extract and implement CPU core APIs in `blst/fft/`.
2. Rewire benchmarks to core APIs.
3. Land correctness tests with oracle-based validation.
4. Add usage documentation and CUDA TODOs.

## Review Checklist
- API shape supports both explicit algorithm selection and precomputed B flow.
- No benchmark-only logic leaked into production API.
- Tests are deterministic and cover size transitions.
- CPU-first scope and CUDA follow-up are explicit.
