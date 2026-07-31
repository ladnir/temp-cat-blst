# CUDA Experimental Policy

## State

CUDA is retained as an experimental, opt-in subsystem. It is not a removal
candidate and is not a primary CatBLST production target. The supported
production profiles remain scalar and AVX-512.

## Policy

- Keep CUDA source, headers, tests, and build integration in the repository.
- Do not claim production readiness, API/ABI stability, or release support for
  CUDA from the scalar/AVX-512 audit results.
- CUDA must remain opt-in and must not be selected silently by supported
  scalar, AVX-512, C, native C++, or Rust configurations.
- CUDA-specific changes, audits, or validation require an explicit CUDA task
  and must record the toolkit, compiler, GPU, driver, and fallback behavior.
- CUDA failures do not block scalar/AVX-512 production readiness unless they
  reveal a defect in shared retained code.
- CUDA benchmarks remain inventoried but excluded and must not be run without
  an explicit benchmark request.

Manifest rows use disposition `out_of_scope_experimental`, status
`experimental`, subsystem `cuda`, and profile `cuda_experimental`.
