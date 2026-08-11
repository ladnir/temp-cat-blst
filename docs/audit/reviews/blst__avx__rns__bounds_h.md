# `blst/avx/rns/bounds.h` Review

## Scope

This record covers the complete compile-time interval algebra through commit
`0471f82`.

## Contracts

- Every interval is closed: `Bounds<L, U>` represents `[L, U]` and requires
  `L <= U`.
- Addition, subtraction, and multiplication compute in the exact signed
  128-bit compatibility carrier before forming signed 64-bit template bounds.
- `recast<L, U>` may only widen an existing interval.
- `shift_left<N>` accepts `0 <= N < 63`, multiplies both endpoints by `2^N`
  in the wide carrier, and rejects a result outside signed 64-bit bounds.
- Bounds are compile-time proof metadata. They perform no runtime validation
  and generate no arithmetic instructions by themselves.

## Findings

- `AUD-239`: the former native-MSVC signed carrier was only 64 bits.
- `AUD-246`: the unused unary negation and the former signed shift expression
  could overflow before producing the intended bound diagnostic.

## Validation

- Compile-time regressions cross both signed 64-bit limits, multiply signed
  boundary values, and verify `[-2, 3] << 2 == [-8, 12]`.
- The serial portable-AVX Debug production/test target compiles every selected
  bound instantiation with the hardened shift contract.
- The MatrixNoK structural validator requires the wide multiplier and forbids
  restoration of the unsafe shift spelling or unused unary negation.
- The focused MatrixNoK arithmetic suite passes 24/24; the combined primitive
  and C API filter passed 32/32 before the final test-only assertion addition.

## Residual risk

The interval algebra proves only its type-level arithmetic. Correct attachment
of an interval to physical AVX lanes is the responsibility of the producing
operation and its higher-level proof or differential tests.
