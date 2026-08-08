# MatrixNoK Construction and Proof Map

## Status and scope

This document is the readable specification of the selected BLS12-381
MatrixNoK candidate. The implementation is performance-selected and is now
semi-frozen for Lean correspondence work, but it is not part of public AVX
dispatch. A proof claim applies only to the exact source and generated-constant
hashes recorded with that claim.

The selected construction is:

1. `matrixnok::FpRingCandidate` for Fp;
2. `FP2Ring<FpRingCandidate>` for Fp2 and G2;
3. `FP12Ring<FpRingCandidate, 3, Fp12ConstantsCandidate>` for Fp12/GT;
4. `G1CurveCandidate` and `G2CurveCandidate` for curve arithmetic; and
5. `MillerCandidate` plus `PairingInverterCandidate` for pairing.

`matrixnok::Fp2RingCandidate` is a retained comparison microkernel. It is not
in this selected call graph, is not dispatched, and is outside the first proof
scope.

## Base-field representation

Each `FpRingCandidate::Element` stores one Fp value in two independent RNS
bases, M1 and M2. Each base has eight pairwise-coprime odd moduli below 2^50;
one AVX-512 lane holds one residue. A `StandardResidue` lane is in `[0, 2m_i)`.
Both bases must decode to the same value modulo the BLS12-381 prime.

`BoundedElement` bounds lane representatives. `ExpandedElement::rns_bounds`
tracks the redundant multiple of the field modulus represented across the two
bases. These types make arithmetic obligations visible to C++ compilation,
but do not validate values at runtime. Lean must prove that every producer
satisfies the bounds claimed by its result type.

Canonical import follows one path:

1. split the canonical integer into eight little-endian radix-2^50 digits;
2. apply the fixed `ConvertToMatrix` into M2;
3. Montgomery-reduce the M2 result; and
4. apply the fixed R2 matrix and M1 reduction to derive M1 from that M2 value.

Export applies `ConvertFromMatrix` to reduced M2 and reconstructs the scalar
limbs. It is an observation boundary, not part of the multiplication kernel.

## Selected multiplication

For standard inputs, the scalar and fixed-batch kernels implement the same five
stages:

1. use IFMA to form products independently in M1 and M2;
2. Montgomery-reduce the M1 product to a standard residue;
3. add correction-free `R1(M1)` to the existing, unreduced M2 product;
4. Montgomery-reduce that combined M2 accumulator; and
5. expand reduced M2 through R2 to reconstruct standard M1.

The fixed bound ledger is:

| Quantity | Bound |
| --- | ---: |
| standard input residue | 2 |
| one elementwise product | 4 |
| eight-term correction-free R1 contribution | 16 |
| conservative delayed-path C++ type allowance | 32 |
| combined M2 product plus R1 contribution | 20 |
| generated base redundant multiple | 40 |
| delayed-addition input accepted by `batch_reduce` | 800 |
| signed redundant-range offset | 1932 |
| maximum tracked redundant multiple | 4309 |

The bound 16 is exactly eight input lanes, each below two source moduli. It is
stronger than the generic delayed-reduction path's retained C++ type allowance
32. The latter is deliberately conservative: narrowing it changes selected
higher-tower code generation and is deferred until a separate benchmark gate.
It is not the mathematical bound Lean should prove. The generated constants
explicitly declare `r1_skip_correction`; the C++ owner statically requires that
declaration. P2 must prove that omitting the quotient correction preserves the
represented value; that the direct multiply's combined bound 20 is below the
generated base redundant multiple 40; and that no 64-bit IFMA accumulator
wraps.

## Fixed matrix semantics

`FixedRNSMatrix` consumes eight source residues. For each source lane it adds a
fixed row times that scalar to independent IFMA high and low accumulators. It
also computes the fixed-point quotient numerator, takes its high 64 bits, and
adds the generated correction vector. Shifted entry points express the same
correction in IFMA radix 2^52.

The high and low output arguments are in/out accumulators:

- conversion callers initialize them to zero; and
- fused MatrixNoK reduction initializes them with an existing product.

This accumulator contract is essential to the correspondence theorem. A model
that treats the arguments as write-only outputs proves a different algorithm.

`r1_change_base_raw` uses the fixed cyclic-diagonal form of the 8x8 R1 matrix.
At diagonal `d`, the source vector has been rotated by `d` lanes and is
multiplied lane-wise by row `d`. The eight products reconstruct the ordinary
matrix-vector product without horizontal reductions. Its input is scratch and
is rotated in place; its high and low outputs are in/out accumulators.

## Field tower layout

The selected Fp2 is `c0 + u*c1` with `u^2 = -1`. Fp12 uses the existing
CatBLST six-Fp2 tower and stores its twelve Fp coordinates in a flat array. Flat
coordinate `i` corresponds to scalar storage
`Fp12[i / 6][(i % 6) / 2][i % 2]`. This ordering is part of the conversion and
differential-test contract; it must not be inferred from memory layout.

The MatrixNoK Fp12 constants policy binds Frobenius constants and delayed
reduction offsets to this exact RNS basis. It contains no runtime state. The
split cyclotomic-square schedule is the selected pairing implementation and is
part of the eventual P3/P4 correspondence.

## Deliberate source duplication

The following duplication is retained because the shapes serve different hot
paths and have separately inspected code generation:

- register-latency and fixed-batch R1 cyclic schedules;
- scalar and fixed-batch base-field multiplication pipelines; and
- the incumbent monolithic and selected MatrixNoK split cyclotomic square.

Do not replace these with runtime abstraction, indirect calls, dynamic loops,
or temporary containers merely to shorten the proof. Lean should define one
mathematical operation and prove each concrete schedule refines it.

## Lean boundary and freeze rule

P0 models the exact lane operations used here. Native intrinsic semantics are
axioms at that boundary; the portable `_cc512_*` implementation and native
differential tests are executable evidence, not the proof itself. P1 proves
the generated constants, dimensions, coprimality, inverses, and bounds. P2
proves the import, multiplication, reductions, and export pipeline. The tower,
curve, and pairing layers follow only after that base-field slice is stable.

After the proof-facing source hash is recorded, changes to a mapped operation,
constant, loop schedule, permutation, bound, or representation reopen its
correspondence record. A reviewed comment-only change may update the source
record without reopening an unaffected theorem.
