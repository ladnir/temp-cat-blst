# CatBLST AVX Arithmetic Proof

This is the machine-checked proof project for the selected BLS12-381
MatrixNoK arithmetic. It implements the layers and claim boundary in
`docs/avx/PROOF_PLAN.md`; `docs/avx/MATRIXNOK_CONSTRUCTION.md` is the readable
specification of the exact selected construction and bound ledger.

The project pins Lean and Mathlib to `v4.32.1`. Updating either pin is a
reviewed dependency change and does not occur implicitly.

## Local setup

Install Elan, then run these commands from `proof/avx`:

```text
lake update
lake exe cache get
lake build
```

Elan selects the exact version in `lean-toolchain`; do not install or select a
different system Lean. `lake-manifest.json` locks Mathlib's transitive
dependencies, while the downloaded `.lake/` cache remains local and ignored.

## First milestone

The first proof slice is intentionally narrow and follows the implementation
dependency order:

1. define 64-bit word, 52-bit digit, IFMA-low, IFMA-high, accumulation, and
   eight-lane permutation semantics;
2. prove the low/high product split and the no-unintended-wrap bounds used by
   the fixed r1 kernel;
3. import the exact generated BLS12-381 constants through a deterministic
   source generator and prove their defining modular relations;
4. prove the fixed cyclic r1 schedule equals the checked-in matrix map;
5. prove both RNS Montgomery reductions and conversion round trips; and
6. prove `FpRingCandidate::multiply` represents scalar multiplication modulo
   the BLS12-381 prime and returns both residue halves below `2m_i`.

No theorem may be marked complete by an executable test, `native_decide`, or a
generated equality alone when the corresponding general arithmetic or bound
claim remains assumed. Computation can discharge fixed finite facts; reusable
arithmetic properties require explicit theorems.

## Correspondence rule

`CORRESPONDENCE.csv` is the review boundary between Lean and C++. Each accepted
row records the exact implementation commit, source hash, constant hash, Lean
declaration, theorem, and native/portable evidence. A mapped C++ or constant
change reopens the row even if tests still pass.

The proof does not establish compiler correctness, x86 ISA correctness, or
microarchitectural constant-time behavior. Those assumptions are tracked by
the native code-generation and constant-time evidence described in the main
proof plan.

## Resource policy

Toolchain installation and Mathlib cache retrieval are explicit setup steps.
Proof builds run serially unless separately approved. They do not build
CatBLST, CUDA, benchmarks, BTx/BTE, or the portable-AVX Release product.
