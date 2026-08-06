# File Review: Pairing context and aggregation implementation

## Record

| Field | Value |
| --- | --- |
| Path | `blst/sign/aggregate.cpp` |
| Roadmap phase / queue | Phase 15 curves and protocols, tenth packet |
| Reviewed commit | `091dfb2cae096ee8ed801036e98827b3253b2e84` |
| Reviewer | Codex audit team |
| Date | 2026-08-05 |
| Disposition | in-scope production |
| Subsystem | batched pairing state machine and signature aggregation |
| Profiles | scalar and selected AVX implementation; C/C++ and Rust boundaries |
| Upstream reference | pinned BLST `src/aggregate.c`, `src/server.c`, and `src/pairing.c` |
| Related files | `blst/sign/aggregate.h`, `blst/c_api.cpp`, `bindings/blst_scalar.h`, `bindings/rust/src/pairing.rs`, `test/Aggregate_audit_tests.cpp` |

## Purpose and trust boundary

This file owns the opaque pairing-context state machine, batches affine pairing
terms, accumulates optional signatures, merges committed worker contexts, and
performs final verification. Messages, augmentations, DSTs, encoded signatures,
public keys, and raw expert points cross trust boundaries here. Public aggregate
operations must report malformed inputs without partially publishing context
state; raw aggregation deliberately requires trusted subgroup affine inputs.
The underlying formulas may evaluate on a wider curve domain, but CatBLST does
not promise that advanced behavior under `AUD-160`.

## Contracts and lifecycle

- Initialization records hash-versus-encode mode and borrows the DST. A DST
  placed immediately after the context uses a control bit and relative recovery,
  so moving an owning allocation does not leave a stale pointer.
- A context records one PK/signature orientation. Calls in the other orientation
  return `BLST_AGGR_TYPE_MISMATCH` without mutation.
- Default aggregation checks a non-infinite signature subgroup, rejects an
  infinite PK, and requires the caller to have subgroup-validated the PK. The
  `chk_n` entry points explicitly select both subgroup checks.
- All scalar-shape, orientation, signature, infinity, and requested PK checks
  complete before orientation, signature, batch, or Miller state is published.
- Pending batches contain at most eight terms. Commit is idempotent when the
  batch is empty. Merge requires both inputs committed and preserves the source.
- The Miller-result getter commits pending terms, returns multiplicative identity
  for an empty context, and exports a pre-final value rather than a GT element.
- Encoded infinity adds the identity: it preserves a supplied accumulator and
  produces identity only when no prior accumulator exists. Exact aliasing and
  distinct input/output objects are supported.
- The implementation allocates no memory. No function in this file throws across
  the C ABI. Weighted aggregation uses variable-time public-scalar operations.

## Compatibility decisions

Ordinary upstream pairing, aggregate, commit, merge, final-verify, and signature
aggregation signatures remain compatible. CatBLST intentionally replaces the
low-level upstream/previous accumulator accessors with
`blst_pairing_get_miller_result(out, ctx)`: caller-owned output avoids a large
by-value C return and avoids adding selected-AVX cache storage to every context.
The semantic name also prevents treating a Miller result as an arbitrary Fp12
or already-finalized GT element.

## Evidence

1. Scalar MSVC Debug and Release each pass six focused tests for empty export,
   trailing DST, G1/G2 infinity accumulation, and byte-exact rejected-item state.
2. Selected portable-AVX Debug builds and passes the same 6/6 tests.
3. Thirteen safe-Rust pairing tests cover owned/empty/growing DST storage, moving
   between threads, merge, scalar-shape rejection, typed orientation, explicit
   omitted signatures, uniqueness handling, and empty Miller identity.
4. Optimized MSVC symbols contain the new C and native getter and no old
   `blst_pairing_as_fp12` accessor.
5. Structural checks find no `dyn Any`, invented pointer 42, aggregate `N_MAX`,
   or old accessor spelling. Public documentation, group-contract, and audit
   inventory checks pass.

## Findings and residual evidence

`AUD-152` through `AUD-159` are validated. No critical or high defect remains
open in this file. Sanitizer/fuzz execution, native AVX-512 hardware execution,
Apple execution, and the complete supported-language milestone matrix remain
release evidence. CUDA and frozen BTx/BTE were neither built nor reviewed. No
benchmark ran because the state-machine schedule and hot pairing kernels did not
change.

## Completion

Complete Phase 15 file review at implementation commit
`091dfb2cae096ee8ed801036e98827b3253b2e84`.
