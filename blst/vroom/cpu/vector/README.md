# cpu/vector

- **AVX512ifma** — wraps AVX intrinsics for arbitrary-length SIMD vectors.
- **`fallback.hpp`** — x86 implementation of the same intrinsics (correctness only, not performance; simulates 52-bit mul with 64-bit ops, shifts, masks).
- **`changebase.hpp`** — core CRNS code; dominates compute time across applications.
- **`multiplication.hpp`** — multiplication via CRNS.
- **`conversion.hpp`** — CRNS conversion to/from radix form.
- **`vector_impl.hpp`** — macros to select AVX512 vs fallback.
