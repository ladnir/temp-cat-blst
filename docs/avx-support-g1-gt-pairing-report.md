# BLST AVX Support Report (G1, GT, Pairing, Batch `n`)

## Scope

This report focuses on AVX support in the current `blst` repo for:

- G1 group addition and scalar multiplication.
- GT group addition and scalar multiplication (additive notation over multiplicative `Fp12`).
- Pairing (`miller_loop` + final exponentiation).
- Batch settings where `n` independent operations are processed together.
- Relevance to FFT-in-the-exponent workflows:
  - forward FFT on `G1` and `G2`
  - inverse FFT on `GT`

## High-Level Summary

- AVX512 infrastructure exists and compiles selected AVX files when `BLST_AVX512` or `BLST_AVX512_PORTABLE` is enabled (`CMakeLists.txt`, `blst/CMakeLists.txt`).
- In-tree AVX group work is concentrated in:
  - `G2` pairing line kernels (`blst/curve/G2_avx.*`)
  - `Fp12` arithmetic kernels (`blst/fields/Fp12_avx.*`)
  - AVX pairing helper implementation (`blst/curve/pairing_avx.*`)
- **Main API paths do not currently dispatch to AVX pairing/group code**:
  - `GT::miller_loop_n` and `GT::final_mult` use scalar/generic paths (`blst/curve/GT.cpp`, `blst/curve/GT.h`).
  - aggregate/signing flows call `GT::miller_loop_n`, not `blst::avx::*` (`blst/sign/aggregate.cpp`).
- `G1` has no dedicated AVX group-op implementation file in `blst/curve`.
- Batch MSM support exists (`multi_scalar*` / `multi_scalar_par*`) but is mostly algorithmic batching; AVX512 SIMD dispatch for G1/GT group ops is not present.

## What Is Present vs Missing

### 1) G1 operations

Present:

- Scalar/generic G1 add/double/mult in `blst/curve/point.cpp`.
- Scalar multiplication uses GLV-specialized path for 256-bit scalars (`G1Point::mult_glsv`) plus generic windows (`Point<Field>::mult`).

Missing for AVX:

- No `G1_avx.cpp` / `G1_avx.h` in `blst/curve`.
- No AVX-dispatched G1 add/double/scalar-mul entry points.
- No G1 AVX tests comparable to `G2_avx_tests`.

Implication:

- Forward FFT on G1 currently uses scalar group ops (`CoeffCtx::plus`, `CoeffCtx::mult` in `blst/fft/context.h` + `blst/fft/ctFft.h`), not AVX group kernels.

### 2) GT operations

Present:

- GT API is implemented via `Fp12` operations (`blst/curve/GT.h`):
  - `GT::add` -> `Fp12::mul`
  - `GT::doubling` -> `Fp12::cyclotomic_sqr`
  - `GT::mult` -> `Fp12::cyclotomic_exp_glv4` (256-bit) or `Fp12::cyclotomic_exp`

Missing for AVX integration:

- GT methods do not dispatch to `Fp12_avx` wrappers.
- AVX `Fp12` functions exist (`mul_fp12_vector`, `sqr_fp12_vector`, `cyclotomic_sqr_fp12_vector`, `mul_by_xy00z0_fp12_vector`) but are used only in `pairing_avx.cpp`.

Implication:

- Inverse FFT on GT currently uses scalar GT operations for butterfly and scaling steps.

### 3) Pairing

Present:

- Mainline scalar pairing path:
  - `GT::miller_loop`, `GT::miller_loop_n`, `GT::final_mult`.
  - `MILLER_LOOP_N_MAX` is 16 in array-based `miller_loop_n`.
  - Pointer-array `miller_loop_n` variant chunks across groups of 16.
- AVX pairing module exists under `namespace blst::avx`:
  - `miller_loop_n`, `final_exp`, `optimal_ate_pairing` (`blst/curve/pairing_avx.*`).
  - `kMillerLoopNMax = 16`.
  - Uses vector-specialized inner kernel only for `n == 1` (`add_n_dbl_n_vector`).
  - Falls back to generic loop for `n > 1` (`add_n_dbl_n_generic`).

Missing:

- No wiring from public/main GT pairing path to `blst::avx::*`.
- No references to `pairing_avx.h` outside its own `.cpp`.
- No dedicated tests that exercise `blst::avx::miller_loop_n/final_exp` directly.

Implication:

- AVX pairing code is present but currently a side path, not the default runtime path.
- For batch `n` pairings, AVX code does not currently provide broad SIMD scaling beyond the `n==1` fast inner path.

## Batch / Parallel (`n`) Assessment

### Pairing batch (`n`)

- Scalar path supports batched pairings:
  - Direct array API limited to `n <= 16`.
  - Pointer-array API handles larger totals by chunking and multiplying chunk results.
- AVX pairing helper also accepts `n` arrays, but the strongest vectorized inner step is currently for `n == 1`.

Net: batch pairing exists functionally, but AVX acceleration for general `n > 1` is incomplete.

### Batch scalar multiplication / MSM

- MSM code exists for G1/G2 (`multi_scalar.*`) and batched pipeline templates for G1/G2/GT (`multi_scalar_par_*`).
- The batched pipeline processes `M` MSMs of `N` points each (`msm_booth_batched_pipeline`).
- CPU branches in `multi_scalar_par_impl.h` are predominantly explicit loops.
- CUDA branches use Thrust/CUB kernels.

Caveat:

- `multi_scalar_par.h` / `multi_scalar_par_impl.h` currently contain `COMPILE_FAIL` markers in non-CUDA branches, indicating this path is CUDA-oriented or currently incomplete for host-only builds when enabled.

## FFT-in-the-Exponent Context

The structure you want is already present:

- Generic FFT template uses group `plus/minus/mult` through context (`blst/fft/ctFft.h`).
- `CoeffCtx` binds:
  - G1/G2 butterflies and twiddle multiplies to group add/scalar-mul.
  - G1xG2 pointwise multiply to pairing via `GT::miller_loop` + `GT::final_mult` (`blst/fft/context.h`).
- `middle_product` code performs forward transforms on G1/G2 and inverse on GT (`blst/fft/middle_product.h`).

What is missing for AVX-first execution:

- G1 forward FFT kernels are scalar (no G1 AVX group kernels).
- GT inverse FFT kernels are scalar at GT API level (no GT AVX dispatch).
- Pairing calls in FFT path use scalar GT pairing path, not `blst::avx`.

## Practical “Coding Effort” Opportunities

### Low to medium effort

- Add compile-time dispatch wrappers in GT pairing entry points:
  - route `GT::miller_loop_n` / final exponentiation to `blst::avx::*` when AVX is enabled.
- Keep scalar fallback unchanged for portability.
- Add tests that compare AVX and scalar outputs for:
  - single pairing
  - batched pairing (`n > 1`)
  - edge cases (infinity inputs)

### Medium effort

- Integrate AVX `Fp12` primitives into GT add/double/final-exp hot paths.
- Add a dedicated GT-AVX benchmark and correctness test matrix.

### Higher effort

- Implement true SIMD-friendly batch pairing for `n > 1` in AVX path (reduce generic fallback usage).
- Add dedicated G1 AVX group-operation kernels (add/double/scalar-mul), then wire FFT contexts to exploit them.

## Bottom Line

For your target workflow (forward FFT on G1/G2, inverse FFT on GT):

- **Present now**: correct algorithms and batching structures, plus AVX kernels for parts of pairing/Fp12/G2 line work.
- **Missing now**: end-to-end AVX dispatch for G1/GT group ops and general-`n` AVX pairing acceleration in the main execution path.
- **Most realistic near-term win**: integrate `blst::avx` pairing/final-exp dispatch and validate in FFT middle-product paths; this should provide benefit without requiring full API redesign.
