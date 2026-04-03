# cpu/vector

- **`avx512ifma.hpp`** — wraps AVX intrinsics for arbitrary-length SIMD vectors through the intrinsic compatibility layer.
- **`changebase.hpp`** — core CRNS code; dominates compute time across applications.
- **`multiplication.hpp`** — multiplication via CRNS.
- **`conversion.hpp`** — CRNS conversion to/from radix form.
- **`vector_impl.hpp`** — unified include for the AVX vector implementation.
