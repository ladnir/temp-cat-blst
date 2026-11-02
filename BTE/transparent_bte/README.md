# Transparent BTE experiment (ePrint 2026/2205)

This directory implements the online algebraic construction and the fast batch-decryption algorithm from [*DKG Is All You Need*](https://eprint.iacr.org/2026/2205.pdf). It is a comparison experiment, not a production cryptosystem.

## Implemented scope

- `SetupFromSharingPolynomial` converts caller-supplied Shamir coefficients into the public and secret values consumed by the online construction. It is a deterministic test and benchmark adapter, not a DKG.
- `Encrypt` produces the four algebraic ciphertext components from Figure 1.
- `PartialDecrypt`, `VerifyPartialDecryptionShare`, and `CombinePartialDecryptionShares` implement the threshold-share path.
- `PrecomputeBatchDecryptionNaive` is the direct quadratic Cauchy transform and serves as a correctness oracle.
- `PrecomputeBatchDecryptionMsm` evaluates each direct Cauchy row with CatBLST's native public-scalar G1 and G2 MSM kernels. In an AVX build, it keeps projective points and Montgomery-form `Fr` values in their native representations and reuses one aligned scratch allocation across all rows.
- `PrecomputeBatchDecryptionQuadraticPairing` is a measured alternative that replaces the quadratic G2 scalar multiplications with prepared pairing lines. It is retained as negative experimental evidence: the wider Miller loops lose to the MSM path except at the smallest measured batch.
- `PrecomputeBatchDecryptionFast` implements Lemma 1 with balanced subproduct trees, group-valued numerator polynomials, scalar and group remainder trees, batch inversion, and two-pair combined Miller loops.
- `FinalizeBatchDecryption` consumes the ciphertext-only precomputation and the combined share. It performs the remaining one pairing per slot.

The implementation uses the selected CatBLST point and field types directly. An AVX-512 build therefore remains in the native AVX representation throughout the computation. Scalar builds use scalar types. Canonical encoding occurs only when `CiphertextPoint` hashes the compressed first component.

## Deliberate exclusions

The paper specifies the relation for a simulation-extractable NIZK but does not instantiate the proof system, transcript, or domain separation needed for an interoperable implementation. The experiment omits that proof from both constructions rather than inventing a favorable implementation. Consequently, this code does not by itself provide the paper's CCA security claim.

The experiment also excludes DKG, wire formats, malformed-point admission, networking, persistence, and production scheduling. Inputs are native, already-validated group elements.

## Algorithm selection

The public experiment exposes all decryption-precompute algorithms explicitly. The crossover is architecture-dependent and therefore is not embedded in the library API. With full-width public scalars on Peach, the native row-MSM implementation beats both repeated multiplication and the Lemma 1 implementation at every measured size from `B=4` through `B=512`. The benchmark harness uses the constant one-slot special case at `B=1`, direct multiplication for the single pair at `B=2`, native MSM from `B=3` through `B=512`, and Lemma 1 for larger, unmeasured sizes. Scalar builds currently use the direct oracle behind the MSM entry point; the optimized comparison reported here is native AVX-512.

Tests compare the direct, MSM, prepared-pairing, and Lemma 1 outputs at irregular sizes and run complete threshold round trips. See [the benchmark report](../../benchmark/TRANSPARENT_BTE_2026_2205.md) for measurements and limitations.
