# `blst/avx/rns/bounded_element.h` Review

## Scope

This record covers `BoundedElement`, `WideElement`, delayed products, and their
complete retained operation surface through commit `0471f82`.

## Contracts

- `BoundedElement<Bounds<L,U>, limbs, element_bits>` uses at least one lane,
  accepts `1 <= element_bits <= 52`, and stores ordinary lane values in the
  signed 64-bit envelope implied by `[L,U]`.
- IFMA multiplication accepts only nonnegative operands whose upper bounds fit
  one 52-bit multiplication digit.
- `WideElement` represents a redundant high/low radix-`2^52` accumulator. Its
  bounds are nonnegative, below 2048, and nondegenerate except for `[0,0]`.
- Wide subtraction is borrow-free: the left lower bound must dominate the
  complete right interval. The high and low vectors are then safely subtracted
  independently without a propagated borrow.
- Default initialization is intentionally overwrite-before-read. `zero()` and
  `WideZero()` explicitly construct zero vectors.
- Recasts only widen metadata. `complete()` materializes delayed IFMA products;
  it allocates no storage and introduces no runtime dispatch.
- Selection, equality, and zero testing inherit the reviewed fixed-control
  vector operations. All methods are internal and forced inline where hot.

## Findings

- `AUD-244`: the generic C++ wide carrier did not independently enforce the
  physical 52-bit accumulator limit already proved for selected MatrixNoK
  paths. Current production's largest bound is 832 and was never incorrect.
- `AUD-245`: borrow-free wide subtraction was documented and proved for current
  formulas but not enforced at its C++ operation boundary.
- `AUD-247`: unused diagnostic conversion/checking helpers contained unsafe
  signed reconstruction and unchecked diagnostic multiplication.
- `AUD-248`: valid radix and initialization contracts were incomplete, and
  zero construction relied on subtle value-initialization behavior.

## Validation

- Every current production and test expression compiles with the new capacity,
  nonnegative-range, nondegenerate-range, and no-borrow static assertions.
- A direct compile-time regression instantiates the largest selected bound
  `832`, pins the capacity to `2048`, and checks `[80,81] - [0,40]` produces
  `[40,81]` under the no-borrow rule.
- A runtime regression confirms `BoundedElement::zero()` produces eight zero
  lanes.
- The serial portable-AVX Debug target builds and the focused MatrixNoK suite
  passes 24/24. The primitive plus C API filter passed 32/32 before the final
  test-only assertion addition.
- The structural validator requires every retained invariant and forbids the
  removed diagnostics and heavy standard-library includes.
- An optimized native MSVC production-Fp probe contains a 312-byte frame,
  eight scalar `mul` instructions, zero calls, and zero conditional branches.
  The current count is 243 instructions versus the prior 244-instruction
  record; no benchmark was warranted.

## Residual risk

The bounds are compile-time claims, not runtime validation of every AVX lane.
Correctness still depends on producers satisfying their declared intervals,
the accepted MatrixNoK proofs and differential tests, and compiler/IFMA
lowering at the existing P0 boundary. Default-constructed storage must be
fully overwritten before read.
