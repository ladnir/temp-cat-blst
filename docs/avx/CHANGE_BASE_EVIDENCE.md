# MatrixNoK Change-of-Base Evidence

## Packet boundary

- Implementation commit: `d8e5adceedd16b81917439d1981d746a81f5bc2e`
- Upstream algorithm commit: `d6ecd534c54a41702ea80ae75b52250b6d51b372`
- C++ owner: `blst/avx/rns/matrixnok_change_base.h`
- Constant owner: `blst/avx/generated_matrixnok_constants.hpp`

This packet adds an isolated BLS12-381 r1 kernel. It has no production caller
and does not change the incumbent AVX field representation or dispatch.

The kernel is specialized to one eight-lane table. Template recursion expands
all eight diagonals, and an index-sequence fold expands every fixed batch item.
It performs seven cumulative cyclic permutations and eight pairs of IFMA high
and low accumulation. It intentionally applies no quotient correction; that
MatrixNoK precondition remains a proof and integration gate.

## Independent executable reference

`AvxMatrixNoKPrimitives.FixedR1PermutationMatchesScalarMatrix` does not call a
second vector implementation. It independently splits each product of two
52-bit values using 26-bit scalar halves, accumulates the low and high IFMA
parts separately with unsigned wrap semantics, and compares every output lane.

Cases cover:

- every input lane at its corresponding modulus minus one;
- small lane-index values with near-`UINT64_MAX` accumulators, exercising
  accumulator wrap;
- deterministic pseudorandom valid residues with nonzero high and low
  accumulators; and
- a fixed batch of three, so source-level batch expansion is exercised.

## Local portable Debug

MSVC 19.50.35728 configured a portable-only Debug profile with native AVX-512
disabled. Only the existing `blst_avx512_compat_tests` target was built, with
one job. The target now supports portable-only execution without adding
another executable.

Focused result: 4/4 passed:

- three `AvxMatrixNoKPrimitives` tests; and
- the hostile fixed-constant validation test.

The existing native-plus-portable MSVC target also compiled the packet. It was
not executed locally because this workstation lacks AVX-512.

## Peach native AVX-512IFMA

Host: AMD Ryzen 9 7950X. Toolchain: GCC 15.2.0, Release, `-march=native`.
Both native and portable support were enabled; CUDA, assembly, and benchmarks
were disabled. Only the existing compatibility target was built, with one job.

The same four focused tests passed natively. Existing `AUD-102`
`-Wignored-attributes` warnings remained visible.

## Optimized code generation

`audit_avx_matrixnok_r1_change_base` was compiled with GCC 15.2.0 using
`-O3 -DNDEBUG -march=native`.

The resulting function contains:

- three input ZMM loads and two final ZMM stores;
- eight fixed `vpmadd52huq` and eight fixed `vpmadd52luq` instructions;
- seven fixed `vpermq` instructions;
- fixed RIP-relative table loads; and
- no branch, call, stack adjustment, conversion object, or spill.

## Residual proof obligations

- Prove that the permuted diagonal schedule equals multiplication by the exact
  checked-in r1 matrix for all lane values, not only tested vectors.
- Prove the IFMA accumulator range and absence of unwanted 64-bit wrap for all
  production input and accumulator bounds.
- Prove the QR representation invariant that permits omission of the quotient
  correction and establish the candidate ring's redundant-range policy.
- Bind the eventual Lean declarations to the exact implementation and constant
  hashes before this kernel receives a production caller.
