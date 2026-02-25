# AVX C++ Refactor Draft (AVX-Only, Scalar Untouched)

## Scope
- Refactor only AVX code paths (`*_avx.*`, AVX layout/helpers, AVX pairing path).
- Do not modify scalar code behavior or scalar API style.
- Keep allocation-free behavior and performance parity as hard constraints.

## Hard Constraints
- No exceptions.
- No `Result` wrapper types.
- No added runtime overhead.
- No C-style array parameters in function signatures.
- Core math/storage types must be trivially default-constructible and trivially copyable.
- No implicit zero-initialization via constructors/member initializers.
- Use `std::span` only on sufficiently high-level APIs; in hot loops, extract raw pointers once and operate on pointers.
- Use `restrict` only where aliasing assumptions are required and valid.
- In debug builds, validate `restrict` assumptions with `BLST_ASSERT` overlap checks.
- Keep/standardize `GT` naming (do not rename to `CyclotomicFp12`).

## Refactor Intent
This is not an ergonomics facade. This is a type/system/signature cleanup to make AVX code match scalar code style and ownership patterns while preserving low-level control.

## Naming Direction

### Type naming
- Use descriptive packed names that preserve schedule encoding.
- Prefer semantic prefixes over anonymous aliases.
- Use `U` prefix for unreduced/widened forms (consistent with scalar `UFp2` idea).

Example mapping (phase-in via aliases first):
- `fp2_4x2x1w` -> `Fp2Pack4x2x1`
- `fp2x2_4x2x1w` -> `UFp2Pack4x2x1`
- `fp2_2x2x2w` -> `Fp2Pack2x2x2`
- `fp2x2_2x2x2w` -> `UFp2Pack2x2x2`
- `fp2_8x1x1w` -> `Fp2Pack8x1x1`
- `fp2x2_8x1x1w` -> `UFp2Pack8x1x1`

### Function naming
- Keep existing math verb patterns (`add`, `sub`, `mul`, `sqr`, `redc`, etc.).
- Keep schedule suffix (`4x2x1`, `2x2x2`, etc.) where needed for clarity.
- PascalCase for types; function naming can remain mixed where practical.

## Type Rules
For each new/renamed AVX type:
- POD only.
- No constructors.
- No default field initialization.
- Add compile-time checks:

```cpp
static_assert(std::is_trivially_default_constructible_v<Fp2Pack4x2x1>);
static_assert(std::is_trivially_copyable_v<Fp2Pack4x2x1>);
```

## API Shape Rules

### High-level entry points
- Accept `std::span` where a sequence is naturally represented (e.g. batch pairing APIs).
- Remove redundant explicit `size` if directly implied by span length.

Pattern:

```cpp
void miller_loop_n(GT& out, std::span<const G2AffinePoint> q, std::span<const G1AffinePoint> p)
{
    BLST_ASSERT(q.size() == p.size());
    auto* __restrict q_ptr = q.data();
    auto* __restrict p_ptr = p.data();
    const std::size_t n = q.size();
    // hot loop uses q_ptr/p_ptr, not span indexing
}
```

### Hot kernels
- Pointer/reference signatures only.
- Use `__restrict` when valid and meaningful.
- Keep output-first convention (`out, in...`).

Pattern:

```cpp
void mul_fp12_vec_v3(Fp2Pack2x2x2& out001,
                     Fp2Pack2x2x2& out101,
                     Fp2Pack2x2x2& out2,
                     const Fp2Pack4x2x1* __restrict ab0,
                     const Fp2Pack4x2x1* __restrict ab1,
                     const Fp2Pack4x2x1* __restrict ab2);
```

## `restrict` + Debug Alias Validation

### Policy
- Add `restrict` only when function semantics already require non-overlap.
- If overlap is semantically allowed, do not add `restrict`.

### Debug checks
Use debug-only overlap assertions where `restrict` is used:

```cpp
#if defined(BLST_DEBUG)
inline bool ranges_overlap(const void* a, std::size_t an, const void* b, std::size_t bn) {
    auto ab = reinterpret_cast<const std::byte*>(a);
    auto bb = reinterpret_cast<const std::byte*>(b);
    return (ab < bb + bn) && (bb < ab + an);
}
#define BLST_ASSERT_NO_OVERLAP(a, an, b, bn) BLST_ASSERT(!ranges_overlap((a), (an), (b), (bn)))
#else
#define BLST_ASSERT_NO_OVERLAP(a, an, b, bn) ((void)0)
#endif
```

## Static Member Functions
Use static member functions when ownership/domain is clear.
- Good:
  - `Fp2::mul_pack4x2x1(...)`
  - `GT::final_exp_avx(...)`
- Keep free functions for cross-domain utilities (layout transforms, low-level lane helpers).

## Invariant Checks (Debug Only)
- Use `BLST_ASSERT` only in debug builds.
- Suggested boundary checks:
  - Affine/Jacobian point form assumptions before heavy loops.
  - Packed schedule shape assumptions at pack/unpack boundaries.
  - Non-alias assumptions for restrict kernels.

## Incremental Plan

### Phase 0: Infrastructure (no behavior change)
- Add alias layer for new type names in AVX headers.
- Add POD static_asserts.
- Add debug overlap helpers/macros.

### Phase 1: Fp/Fp2 AVX signature cleanup
- Convert function params from C arrays to explicit type refs/pointers.
- Keep implementation logic identical.
- Introduce `restrict` where already semantically required.

### Phase 2: Fp4/Fp6/Fp12 AVX cleanup
- Apply naming/type migration to packed and unreduced intermediates.
- Preserve algorithm and schedule behavior.

### Phase 3: G2/Pairing AVX cleanup
- Update AVX pairing interfaces to span at high level and pointers in hot loops.
- Keep GT naming and existing behavior.

#### Phase 3 Staging (current)
- Added AVX-only high-level span overloads:
  - `blst::avx::miller_loop_n(GT&, std::span<const G2AffinePoint>, std::span<const G1AffinePoint>)`
  - `blst::avx::optimal_ate_pairing(GT&, std::span<const G2AffinePoint>, std::span<const G1AffinePoint>)`
- Existing pointer+`n` kernels remain the hot-path implementation.
- Added defensive `n==0` handling in AVX `miller_loop_n` entry to avoid undefined behavior.
- Next implementation slice:
  1. Move `blst::avx` pairing APIs from C-array syntax to pointer/reference form where still needed.
  2. Add `restrict` only on AVX hot loops with clear non-alias requirements and debug overlap assertions.
  3. Keep scalar `GT::*` paths untouched; AVX dispatch integration remains additive and explicit.
  4. Validate with:
     - `ctest --test-dir out/build/avx512-portable-check --output-on-failure -R "Avx|avx"`
     - `ctest --test-dir out/build/avx512-portable-check --output-on-failure`

### Phase 4: Remove transitional aliases
- After AVX compile/tests/bench parity, delete old alias names in AVX code.

## Acceptance Criteria Per Phase
- Builds pass for AVX-enabled targets.
- Pairing and field tests pass.
- Benchmarks show no regression (target: parity or better).
- No scalar code changes required.

## Suggested First Patch Set
1. Add typed aliases + POD assertions in AVX headers.
2. Convert one narrow hot path (`Fp2_avx` kernels) to new signatures.
3. Add debug overlap assertions for newly `restrict`-annotated kernels.
4. Verify compile + AVX tests + benchmark subset.

## Notes
- This plan intentionally avoids deprecation machinery because there are no external users yet.
- This plan intentionally avoids higher-level wrappers/facades.
