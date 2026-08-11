# `blst/avx/rns/ring_element.h` Review

## Scope

This record covers the complete retained one-basis and two-basis RNS expression
carrier through commit `1e0b6e1`.

## Contracts

- `SmallElement<Element, RNSBounds>` carries one physical RNS basis.
  `ExpandedElement<Element, RNSBounds>` carries two bases that must encode the
  same field value.
- `RNSBounds` is compile-time proof metadata for the represented multiple of
  the field modulus. Constructors attach that claim but do not validate lanes,
  cross-base consistency, or field semantics at runtime.
- Default initialization is overwrite-before-read. Value initialization with
  `{}` zero-initializes the selected trivial carrier types.
- Addition, subtraction, fixed scaling, shifting, and recasting propagate
  physical and RNS bounds independently through their reviewed lower-level
  operations.
- Multiplication creates two delayed products. `complete()` materializes them
  as wide accumulators. The specialized addition overload is intentionally
  ordered: its left operand is an existing wide accumulator and its right
  operand contributes one delayed product.
- `ctime_select` and `ctime_replace_if` operate only between identical physical
  and proof-metadata types and inherit the fixed-control vector selection.
- `is_zero` and `equal` compare exact RNS storage. They are not semantic field
  predicates when redundant representatives are permitted. Field observation
  boundaries must normalize or convert, as established by `AUD-225`.

## Findings

- `AUD-249`: the unused `from_constant` helper attached `[0,1]` and cross-base
  claims to arbitrary arrays without validation; the unused wide-zero helper
  unnecessarily enlarged the foundational surface.
- `AUD-250`: `complete()` leaked `const` into the proof-metadata template type
  through `decltype(rns_bounds)`.
- `AUD-251`: the small-carrier recast accepted 32-bit bounds while every other
  bound interface uses `int64_t`.
- `AUD-252`: the narrow construction, initialization, delayed-operation, and
  raw-observation contracts were absent from the owning header.

## Validation

- The complete serial portable-AVX Debug target compiles every current
  expression after the stronger and narrower surface changes.
- `RingCarrierContractsAreExact` checks 64-bit RNS recasting, exact non-`const`
  completion metadata, delayed product materialization and accumulation,
  value-initialized zero, exact-representation predicates, and both conditional
  selection directions.
- The structural validator requires the retained contracts and rejects the
  deleted helpers, 32-bit recast signature, and cv-qualified completion form.
- All 20 accepted AVX correspondence rows and 59 mapped hashes still validate.
  This generic carrier is not itself a direct mapped theorem owner.

## Residual risk

Correctness depends on reviewed producers attaching truthful physical and RNS
bounds and constructing mutually consistent bases. Exact predicates remain
deliberately narrower than semantic field observation. Native compiler and
IFMA lowering remain the existing P0 boundary.
