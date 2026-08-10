# cpu/vector

- **`avx512ifma.hpp`** — fixed-register vector operations routed through the selected native or portable intrinsic layer.
- **`vector_impl.hpp`** — unified include for the AVX vector implementation.

The selected BLS12-381 MatrixNoK reduction and conversion schedules live under
`blst/avx/rns`. CatBLST does not retain the earlier generic runtime CRNS engine.
