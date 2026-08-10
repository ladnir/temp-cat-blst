# MatrixNoK Primitive Evidence

## Packet

- Port map: `44fc26f`
- Primitive implementation: `d7e9a71`
- Code-generation probes: `31b647d`, `33a83ee`
- VROOM algorithm source: `d6ecd534`

This packet does not select MatrixNoK or change the incumbent field
representation. It establishes only the BLS12-381-reached vector, conversion,
and bounded-reduction primitives.

## Local portable Debug

Toolchain: MSVC 19.50.35728, C++ latest mode, portable AVX selected.

- `intrin_compat.cpp` and `Avx_matrixnok_primitive_tests.cpp` compile as focused
  objects with one Ninja job.
- The focused temporary link ran
  `AvxMatrixNoKPrimitives.*`: 2/2 passed.
- The combined native-plus-portable compatibility target compiles with the new
  source under the existing unsupported-MSVC probe configuration.

The ordinary AVX-default Debug preset also exposed `AUD-210`, an unrelated
pairing-context static assertion introduced by the earlier lifetime repair.
That failed configuration is not counted as primitive evidence.

## Peach native AVX-512IFMA

Host: AMD Ryzen 9 7950X. Toolchain: GCC 15.2.0, Release, `-march=native`.
Configuration enabled both `BLST_AVX512` and `BLST_AVX512_PORTABLE`, disabled
CUDA, assembly, and benchmarks, and built only
`blst_avx512_compat_tests` with one job.

Result: 12/12 passed:

- 10 direct `_cc` versus `_mm` compatibility tests, including low-128
  extraction, addition, shifts, and broadcast; and
- 2 selected-vector MatrixNoK primitive tests, including reduction gaps up to
  the new 40x bound.

## Optimized code generation

`test/audit_avx_matrixnok_codegen.cpp` was compiled with GCC 15.2.0 using
`-O3 -DNDEBUG -march=native`.

Lane-zero extraction lowers to a ZMM load followed by `vmovq`; it has no stack
round trip, conversion object, wrapper call, or out-of-line helper.

The reduction-quotient broadcast lowers to:

```text
vpsrlq $52
vpaddq
vpbroadcastq
vmovdqu64
```

The bound-40 reducer lowers to six fixed `vpsubq`/`vpminuq` pairs followed by
one store. It contains no conditional branch, call, stack adjustment, or
variable-indexed access.

## Residuals

- `AUD-102` is validated by the later `d222e05` integration packet. Direct
  native-vector `std::array` spellings now use an exact alias whose GCC-only
  diagnostic suppression is confined to the alias declaration. Static
  assertions preserve element size, aggregate size, and alignment; GCC 15.2
  emits direct ZMM code with `-Werror=ignored-attributes`.
- The compile-time reduction coverage theorem is specified for formalization
  in `PROOF_PLAN.md`; the current evidence is exhaustive over selected boundary
  quotients, not yet machine-checked mathematics.
- MatrixNoK constants, permutation arithmetic, RNS correspondence, field
  operations, and performance remain future packets.
