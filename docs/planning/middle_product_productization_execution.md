# Middle Product Productization Execution Notes

## Completed
- Added a reusable CPU API in `blst/fft/middle_product.h`:
  - `MiddleProductAlgorithm` with `BasicFFT` and `SpecializedFFT`.
  - `middle_product_fft_size(...)` helper.
  - `precompute_middle_product_b(...)` for offline FFT(B) preparation.
  - `middle_product(...)` overloads for raw `B` and precomputed `B`.
- Rewired benchmark entry points to the new API for productized calls:
  - `bench_middle_product_basic`
  - `bench_middle_product_fft`
- Added API-level correctness tests in `test/Fft_tests.cpp`:
  - Basic path vs `O(n^2)` field middle-window oracle.
  - Specialized path vs same field oracle.
  - Raw-`B` vs precomputed-`B` equivalence and cross-algorithm equivalence.

## CPU/CUDA Status
- CPU implementation is now productized and directly callable.
- CUDA parity is not implemented yet for the new API.
- A TODO marker is present in `blst/fft/middle_product.h` for GPU parity:
  - `TODO(BLST_CUDA): mirror middle_product() and precompute_middle_product_b() on GPU.`

## Validation Status
- CMake configure succeeds with `x64-Release`.
- Full build/test execution is currently blocked in this environment by missing CRT/Windows SDK headers during compile (`stdio.h`, `corecrt.h`), even when launched through `C:\agent\bin\msvc-env.cmd /vs:17`.
