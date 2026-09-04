# BTX Unfreeze Record

## State

The BTX freeze was lifted by explicit project-owner direction on 2026-08-17.
The current BTX snapshot may be benchmarked, reviewed, modified, or extracted.
This decision does not make the snapshot production-ready.

PFE, TRX, and unrelated shared BTE code remain under their prior classification
until the project owner changes their scope.

## Implementation decision

The first equivalent benchmark comparison completed on 2026-08-17. The result
does not justify porting the complete C++ protocol implementation. Rust owns
BTX-internal protocol types, validation, and cryptographic orchestration. C++
supplies narrow kernels when equivalent measurements show a material
advantage. Consensus-facing integration is owned elsewhere.

The canonical active plan is
`monad-bft-private/monad-encrypted-tx/docs/BTX_WORK_PLAN.md`. This repository
does not keep a duplicate roadmap.

## Active C++ work

- Admission and combination-validation parity is complete.
- G1 MSM arithmetic parity is implemented and passes the shared vectors. Its
  performance gate remains open because the current laptop measurements are
  too variable for a stable sub-10% comparison at every batch size.
- Batched Miller-loop verification is complete in C++ and implemented for
  review in Rust. The Rust wrapper uses the same `batch_size + 1` equation and
  is 3.0% to 16.2% faster in the final sequential measurements.
- Scalar prepared-alpha comparison is complete. Rust selected the existing
  generic prepared Miller-line boundary rather than ordinary `update`; no
  BTX-specific C++ kernel was added.
- A stateless native-AVX prepared-pairing range kernel is implemented for
  ciphertext-side verification precomputation. It returns a typed pre-final
  Miller product, allocates no memory, and reuses the prepared G2 table already
  retained for alpha processing. Peach evidence is recorded in
  `BTX_PREPARED_PAIRING_RANGES_2026-09-03.md`.
- The matching scalar prepared-G2 range kernel is implemented for review. It
  reuses one immutable scalar line table, batch-converts native projective G1
  ranges, and returns composable pre-final Miller products without evaluation
  allocation. Sixteen-pair serial partitioning has no measurable overhead;
  four-worker scaling reaches 3.82x at 256 pairs. Evidence is recorded in
  `BTX_SCALAR_PREPARED_PAIRING_RANGES_2026-09-03.md`.
- Batch 5 malformed-input parity is complete. The C++
  public-parameter decoder checks the exact SRS dimension and all encoded
  element counts before allocation. Threshold-request decoding likewise
  rejects an impossible element count before reserving storage. Rust owns the
  broader fallible protocol interfaces that have no C++ counterpart.
- Keep the legacy KEM-only ciphertext separate from the full Rust protocol
  type until an approved plan requires convergence.
- Record benchmark profiles, commands, host details, and raw results.
- Never run two benchmark processes concurrently.

Context-separated zero sharing remains blocked on DKG integration. Associated
data, sponsor binding, release evidence, tagged decryption outcomes, and
execution artifacts remain blocked on the primary author's specification
review.

## Production boundary

The C++ and Rust implementations may coexist during kernel comparison. A
passing vector or benchmark establishes only the behavior named by that test.
Neither implementation receives a production-readiness claim until the
specification review and applicable validation gates are complete.
