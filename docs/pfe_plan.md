# PFE Implementation Plan

## Goal

Implement a new `BTE/pfe` module that follows the partial-fraction prior work from
`docs/2026-674.pdf`, but exposes a practical C++ KEM-style API for this repository.

The first milestone is correctness. Performance work comes after the algebra is
implemented and covered by tests.

## Scope For This Pass

- Add a new sibling module at `BTE/pfe/`.
- Focus on C++ only.
- Omit serialization and C bindings for now.
- Omit NIZK / robustness for now.
- Implement threshold share generation and combination in the same style as `btx`.
- Use a fixed batch size `ell` per setup, as in the paper.
- Return per-ciphertext symmetric keys, not encrypted payloads.

## Paper Mapping

The paper defines a batch encryption scheme over partial fractions

- `p_i(X) = 1 / (X + i)`

with setup sampling a secret `x`, publishing:

- a G1 encryption key built from `sum_{i=-1}^ell [p_i(x)]_1`
- public decryption key elements `[p_j(x)]_2` for `j in [ell]`

and threshold-sharing each scalar `p_j(x)` across servers.

The paper's threshold batch encryption construction then computes:

- per-server share:
  `sbk_j = sum_{i in [ell]} s_{i,j} * ct_i[1]`
- combined pre-decryption key:
  `sbk = sum_{i in [ell]} [r_i * p_i(x)]_1`

and decrypts each ciphertext using the `L_{j,1}` / `L_{j,2}` equations from
Construction 1 / Construction 2.

## KEM Adaptation

The paper encrypts a payload by including a third ciphertext component that masks
the message with `r * [p_0(x)]_T`.

For this repository we only need a KEM output. The adaptation for this pass is:

- ciphertext becomes `(c1, c2)` where:
  - `c1 = [r]_1`
  - `c2 = r * sum_{i=-1}^ell [p_i(x)]_1`
- the encapsulated secret is:
  - `R = [r * p_0(x)]_T`
- the exported symmetric key is:
  - `K = HKDF(serialize(R))`

This preserves the core decryption algebra from the paper while matching the
existing `btx`/`trx` usage in this repo, where higher layers handle message
encryption outside the primitive.

To make encapsulation efficient, the public parameters for this implementation
will publish:

- `sum_p_g1 = sum_{i=-1}^ell [p_i(x)]_1`
- `p0_gt = [p_0(x)]_T`
- `dk_j = [p_j(x)]_2` for `j in [ell]`

The `p0_gt` publication is an implementation convenience for the KEM interface.
It matches the paper's note that encryption can avoid a pairing by publishing
the target-group value directly.

## API Shape

The first C++ API should mirror the current `btx` flow:

- `Setup`
- `Encrypt`
- `BuildThresholdRequest`
- `CreatePartialDecryptionShare`
- `CombinePartialDecryptionShares`
- `PrecomputeBatchDecryption`
- `FinalizeBatchDecryption`
- `RecoverBatchKeys`

Planned core types:

- `PfeParameters`
- `PublicParameters`
- `ServerSecretShare`
- `SetupResult`
- `Ciphertext`
- `EncryptionResult`
- `ThresholdRequest`
- `PartialDecryptionShare`
- `BatchDecryptionPrecomputation`

## Fixed-Batch Semantics

The paper is naturally fixed-length in `ell`, with smaller batches handled by
dummy padding. For this implementation:

- every decryption batch must have exactly `ell` ciphertexts
- callers that want smaller logical batches must pad externally

This is faithful to the paper and keeps the implementation simpler.

## Correctness-First Decryption Plan

Even in the correctness-first phase, we can implement the paper's amortized
decryption equations rather than the naive quadratic-pairing variant. The work
split is:

1. Setup-time precomputation in public parameters:
   - `sum_dk = sum_j dk_j`
   - per-slot weighted G2 sums
     `W_j = sum_{i != j} (1 / (j - i)) * dk_i`
   - per-slot self coefficient
     `s_j = sum_{i=-1, i != j}^ell 1 / (j - i)`
   - per-slot coefficient table `1 / (j - i)`

2. Batch precomputation:
   - `tm_j = e(c1_j, dk_j)`
   - weighted ciphertext sums
     `C_j = sum_{i != j} (1 / (j - i)) * c1_i`
   - weighted target-group sums
     `T_j = sum_{i != j} (1 / (j - i)) * tm_i`

3. Finalization:
   - compute `L_{j,1}` and `L_{j,2}`
   - solve `R_j = j * ((j + 1) * L_{j,1} - L_{j,2})`
   - derive the batch key for slot `j` from `R_j`

This keeps us close to the paper and avoids baking in a clearly temporary
quadratic-pairing design.

## NIZK Status

The paper's full CCA-secure construction adds a simulation-extractable NIZK.
The proof relation is essentially:

- witness: randomness `r`
- statement: `(ct[1], ct[2])`
- constraints:
  - `ct[1] = [r]_1`
  - `ct[2] = r * ek[1]`

The paper later instantiates this as a Schnorr-style SE-NIZK with an extra
linear relation. We are intentionally skipping this in the first pass.

## Initial Deliverables

- `BTE/pfe/types.h`
- `BTE/pfe/core.h`
- `BTE/pfe/core.cpp`
- CMake wiring for the new module
- unit tests in `test/Pfe_tests.cpp`
- this planning doc

## Follow-Up Work

- Add serialization support.
- Add C API / bindings.
- Add robust combine / share verification.
- Add the paper's proof system if CCA-hardening is needed in-tree.
- Add benchmarks against `btx`.
