# File Review: scalar Fp addition chains

## Record

| Field | Value |
| --- | --- |
| Paths | `blst/scalar/fields/recip-addchain.h`; `blst/scalar/fields/sqrt-addchain.h` |
| Roadmap phase / queue | Phase 14 scalar fields, generated-addition-chain packet |
| Reviewed implementation commit | `529a79377b6ef848ec90fce35e7580a9f451fc4f` |
| Reviewer | Codex audit team |
| Date | 2026-08-06 |
| Disposition | in-scope production |
| Subsystem | Fixed-exponent BLS12-381 base-field inversion and reciprocal square root |
| Profiles | scalar host implementation; selected/native C++ Fp and Fp2 callers |
| Generation/provenance | `kwantam/addchain` Bos-Coster width-four schedules translated from upstream BLST |
| Upstream references | reciprocal `79400f963d1464c53995d388ef1dcc38edbb8ebf`; reciprocal square root `d02b0d86e25d700ce4f4cb0eccac1f743aea6986` |
| Related files | `blst/scalar/fields/Fp.h`, `blst/scalar/fields/Fp2.h`, `scripts/audit/check_addchains.py`, `test/Addchain_audit_tests.cpp` |

The two generated headers are reviewed exactly at the implementation commit
above. Broader milestone evidence is tracked separately from this file review.

## Purpose and trust boundary

The headers contain fixed, generated exponentiation schedules for p-2 and
(p-3)/4 in BLS12-381 Fp. The first implements fixed-time Fermat inversion for
nonzero values while deliberately mapping zero to zero. The second feeds Fp
square-root, inverse-square-root, and map-to-group operations. Operands may be
secret; the schedule and memory access pattern must not depend on field values.

These are generic arithmetic primitives, not secret containers. By explicit
user decision they do not attempt incomplete transient stack scrubbing. The
library erases narrowly identified keys, scalar decompositions, schedules, and
boundary-owned buffers; comprehensive compiler-created arithmetic spill
erasure is not promised by the constant-time contract.

## Contracts and invariants

- Only scalar BLS12-381 Fp may instantiate either schedule.
- Reciprocal computes `inp^(p-2)`: nonzero input is inverted and zero maps to
  zero for branchless projective-infinity normalization.
- Reciprocal square root computes `inp^((p-3)/4)`, including zero mapping to
  zero. `Fp::recip_sqrt(0)` succeeds with output zero; for nonzero success,
  `out^2 * inp == 1`.
- Exact output/input aliasing is supported because the input is copied into a
  typed live scratch object before output publication. Partial overlap is not
  supported by the owning Fp contract.
- Both schedules have fixed operation and address sequences, allocate nothing,
  perform no runtime dispatch, and invoke only audited Fp square/multiply
  primitives.
- The public C fixed-inverse entry retains its narrower nonzero precondition;
  the wider internal zero behavior exists for branchless point normalization.

## Exact schedule proof

`scripts/audit/check_addchains.py` symbolically interprets every active typed
assignment, square, multiply, and combined square/multiply. It independently
pins the target exponents, operation counts, input slots, and originating BLST
commits. It rejects uninitialized reads, mismatched hexadecimal traces,
unrecognized field operations, multiple output writes, global/generic scope
regression, typed-copy regression, and any final exponent or count change.

The proof results are:

- p-2: 379 squarings, 81 multiplications, 461 atomic trace steps;
- (p-3)/4: 377 squarings, 80 multiplications, 458 atomic trace steps;
- all 16 generated scratch slots are initialized before use;
- every recorded intermediate hexadecimal exponent is exact.

## Checklist evidence

1. Both complete generated files, owning Fp/Fp2 calls, signing, map-to-group,
   public C exposure, and focused tests were traced.
2. Exponent, zero, alias, type-domain, timing, allocation, and output contracts
   are explicit.
3. Original BLST commits and addchain commands are recorded in each header.
4. Every arithmetic step and trace comment is symbolically verified.
5. Generic exponentiation independently checks deterministic boundary values.
6. Exact output aliasing agrees with the independent result.
7. The schedules contain no field-dependent branch or indexed access.
8. No broad scratch erasure was added; the security boundary is documented.
9. Native Fp/Fp2 zero behavior now matches signing and map-to-group callers.
10. Helpers live only in `blst::scalar::detail` and reject non-Fp types.
11. There is no allocation, exception path, threading, or runtime dispatch.
12. Scalar MSVC Debug and Release focused execution passes.
13. Optimized MSVC frame sizes and arithmetic call counts are unchanged; the
    typed assignment removes one entry pointer-comparison branch per chain.
14. Sanitizer, fuzz, Apple, and broad supported-matrix execution remain
    milestone evidence rather than change-loop work.

## Findings

`AUD-181` through `AUD-184`. Proposed `AUD-180` stack erasure was explicitly
declined: generic arithmetic will not pay for an incomplete transient-erasure
claim.

## Validation

- The retained symbolic verifier passes both exact schedules.
- Scalar MSVC Debug and Release each pass 3/3 focused tests covering generic
  exponentiation agreement, zero, one, small values, p-1, exact aliasing, Fp2
  zero, and reciprocal-square-root zero.
- Those tests now belong to the ordinary `blst_tests` source set. The 3/3
  executions used the identical source in a temporary focused target before
  consolidation; scalar and portable-AVX configurations accept the consolidated
  source graph, while a post-consolidation executable run remains deferred.
- Optimized MSVC retains reciprocal/sqrt frames of 880/832 bytes and call
  counts of 86/85. There is no allocator reference or field-dependent branch.
  Typed assignment removes the previous entry pointer-identity branches.
- No arithmetic schedule, benchmark, CUDA, BTx/BTE, portable-AVX Release,
  sanitizer suite, Apple CI, or broad test sweep ran.

## Residual risk and decision

The exact headers are complete at the recorded implementation commit. The
verifier proves the exponent schedule relative to the audited Fp
operations; it is not a formal proof of the assembly multiplication backend.
Constant-time arithmetic does not promise complete erasure of compiler-created
register or stack transients. A stronger erasure threat model would require a
separate end-to-end secret-operation design rather than generic per-function
zeroing.
