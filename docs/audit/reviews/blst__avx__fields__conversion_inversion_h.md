# `blst/avx/fields/conversion_inversion.h` Review

## Scope

This record covers the complete retained MatrixNoK scalar-conversion and
inversion boundary through implementation commit `dd4e4a1`.

## Contracts

- Scalar `Fp` values crossing this boundary use CatBLST's internal reduced
  Montgomery representation. MatrixNoK values use the selected eight-lane,
  two-basis standard RNS representation.
- `ring_to_scalar_fp` consumes a valid selected-Ring value. Raw high/low
  radix-`2^52` accumulators are implementation intermediates and are not field
  values. Export is allocation-free, nonthrowing, and constant-time with
  respect to the represented value, subject to the separately inspected
  scalar assembly kernels.
- `scalar_fp_to_ring` consumes a valid internal Montgomery `Fp` and returns the
  congruent selected-Ring standard element without validation or allocation.
- Fp12 bridges preserve the established flat coefficient order used by GT and
  FFT. They perform twelve independent base-field conversions.
- `invert_fp_vartime_via_scalar` and `ScalarFpVartimeInverter` require a
  nonzero Fp value. `FP12RNSVartimeInverter` requires a nonzero Fp12 value.
  These operations are allocation-free and nonthrowing but variable-time.
  They are suitable for public point normalization and final exponentiation,
  not secret-dependent inversion.
- Secret projective normalization bypasses the variable-time providers and
  calls the fixed-chain scalar reciprocal directly, then erases scalar
  temporaries.

## Findings and remediation

- `AUD-259` repairs secret-dependent branches in the byte-weaving carry chain.
  The established add-with-carry primitive now lowers to `adc`/`setb`, and the
  conversion boundary is forced inline.
- `AUD-260` removes stale mutable references and redundant complete-element
  copies throughout Fp and Fp12 export.
- `AUD-261` removes unused non-Montgomery modes, generic fallbacks, the direct
  scalar-Fp12 inverter, and the unused convenience pairing wrapper.
- `AUD-262` makes variable-time and nonzero contracts explicit in names,
  documentation, and Debug assertions.
- `AUD-263` independently bounds every selected raw-export accumulator,
  adjacent join, and the final REDC input.

The file now carries Category Labs attribution, Apache-2.0 SPDX metadata, and
the VROOM algorithm-lineage note.

## Validation

- `validate_matrixnok_constants.py --self-test` independently models maximal
  residues below `2*m2[i]`. The quotient is 54 bits; raw lanes are at most 55
  bits; joined digits are at most 56 bits; and the reconstructed 435-bit value
  is strictly below the 765-bit `p*2^384` REDC bound.
- `validate_matrixnok_dispatch.py` pins the selected static conversion paths,
  const contracts, explicit variable-time providers, Debug final-carry check,
  and absence of every removed fallback or legacy spelling.
- A serial two-worker portable-AVX Debug rebuild succeeds. Thirty-two focused
  MatrixNoK, GT, inverse-FFT, Miller, and pairing tests pass.
- An optimized MSVC 14.50 NDEBUG AVX-512 probe shows the selected export fully
  inlined. Carry propagation uses `adc` and `setb`; there is no secret-dependent
  branch, address, or conversion-wrapper call. Remaining loops have public
  fixed counts. The 648-byte inlined frame is smaller than the prior nested
  360-byte caller plus 432-byte callee frames.

## Residual risk

The fixed constants and source-level arithmetic are independently bounded, but
x86 IFMA semantics and the scalar assembly kernels remain accepted compiler/ISA
boundaries. Clang native lowering should be rechecked on Peach at the next
authorized native evidence run. No benchmark, remote job, CUDA, portable-AVX
Release build, sanitizer sweep, broad suite, or BTx/BTE audit ran.
