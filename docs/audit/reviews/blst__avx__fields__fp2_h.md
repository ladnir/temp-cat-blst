# `blst/avx/fields/fp2.h` Review

## Scope

This record covers the complete retained selected AVX Fp2 carrier and ring
owner through commit `65206c9`.

## Contracts

- `FP2<X,Y>` stores `x + u*y` with `u^2=-1`. Coordinate template types carry
  the physical and represented-multiple bounds; constructors do not validate
  those producer claims.
- Default initialization is overwrite-before-read. Value initialization
  zeroes the selected trivial coordinate carriers.
- `prep_left` prepares `(x,y,-y)` for reuse. `FP2N` is the prepared-left
  carrier; it is not a distinct field domain.
- Multiplication, squaring, accumulation, fixed scaling, nonresidue
  multiplication, negation, and batch reduction are allocation-free and have
  fixed schedules with respect to field data.
- Variadic batches are nonempty and use the strided layout
  `[x0,...,xN,y0,...,yN]`. Batch size is a compile-time constant.
- `FP2Ring` owns its base-ring policy by value. The selected empty ring occupies
  no additional storage through `no_unique_address`; it cannot dangle.
- `is_zero` and `equal` are exact-representation predicates. Semantic field
  observation must normalize or export redundant RNS values.

## Findings

- `AUD-264`: a second comparison-only Fp2 microkernel duplicated arithmetic,
  tests, probes, and private base-field hooks without a production caller or
  formal correspondence.
- `AUD-265`: `FP2Ring` stored a borrowed `const Ring&` even though construction
  accepted temporaries, permitting a dangling owner.
- `AUD-266`: dead batch plumbing and an implicit prepared-to-unprepared
  conversion enlarged the internal surface; one read-only accumulator input
  was mutable and base-field scaling materialized an avoidable Fp2 temporary.
- `AUD-267`: bound reattachment bypassed the exact metadata alias, provenance
  was missing, and raw predicate and ownership contracts were not stated at
  the Fp2 owner.
- `AUD-278`: the prepared-left carrier lacked a trivial default constructor,
  forcing prepared-line storage to run reduction-bearing zero preparation.

## Validation

- The serial portable-AVX Debug production/test target compiles the selected
  owner. Three focused scalar-oracle, completion-boundary, and fixed-context
  tests pass.
- The structural validator rejects the deleted microkernel, friendship,
  private hooks, borrowed ring, dead helper, implicit conversion, and stale
  proof-type spelling.
- The small MSVC native optimized probe compiles in isolation. Multiplication
  and square contain zero calls; their five branches are fixed-control
  matrix/permutation loops.
- Value-initialized prepared-line tables now compile to one exact zero-fill;
  no per-line reduction, conversion, temporary, or loop remains.
- Accepted `AVX-P3-001` still proves the exact formulas, tight bounds, strided
  layout, and completion pipeline. All 21 correspondence rows and 61 mapped
  hashes validate after exact rebinding.

## Residual risk

x86 IFMA and compiler lowering remain the accepted P0 boundary. Full
sanitizer, fuzzing, statistical leakage, and platform evidence remain milestone
work. No benchmark, remote run, broad suite, CUDA, BTx/BTE audit, or
portable-AVX Release build ran for this packet.
