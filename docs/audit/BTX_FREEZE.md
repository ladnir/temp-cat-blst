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
- Expose batched Miller-loop and prepared `update` operations through narrow,
  allocation-conscious boundaries when the Rust batches reach those stages.
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
