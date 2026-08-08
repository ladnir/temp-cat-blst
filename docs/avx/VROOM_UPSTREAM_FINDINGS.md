# VROOM Upstream Findings

This ledger tracks confirmed defects or hardening changes in code derived from
VROOM that should be offered upstream after CatBLST has a reviewed repair. It
is separate from CatBLST's production audit: every entry links to an `AUD-*`
finding, but this file records the smaller upstream-facing reproducer, affected
VROOM paths, and proposed pull-request shape.

Do not submit a report based only on a suspected correspondence gap. An entry
becomes `ready` only after CatBLST has a concrete reproducer, a minimal repair,
and correctness evidence suitable for upstream review.

## Ledger

| ID | CatBLST finding | Severity | State | Summary | VROOM paths | Upstream action |
| --- | --- | --- | --- | --- | --- | --- |
| VROOM-001 | `AUD-221` | high | fixed in CatBLST; upstream packet pending | MatrixNoK ConvertTo can select a negative CRT representative, violating the immediately following conversion/expansion contracts | `cpu/precompute/precompute.hpp`, `src/bounded_ring.hpp`, generated BLS12-381 MatrixNoK parameters | Submit the canonical witness, focused import/round-trip and multiply-by-one regressions, and CatBLST's floor-estimate repair |

## VROOM-001: ConvertTo may emit a negative representative

### Observable failure

The following valid canonical BLS12-381 field element is below `p`:

```text
1267253245981839142264863935160608057326608172767180188451362181929717619810821363803276830181705476116338914753674
```

Its six little-endian 64-bit limbs are:

```text
11721118646859857034, 14515276896522647911,
16519252582114724462, 14272158255930264238,
5025824079061901589, 593286955712011240
```

Its eight little-endian radix-`2^50` digits are:

```text
500616628141194, 1029601075259562, 686812189624822,
769198589660543, 581522751939018, 26173626664305,
89497429892085, 552541535
```

The generated ConvertTo ceiling estimates choose the exact CRT quotient plus
one. The corrected representative is therefore negative. The M2 Montgomery
reduction preserves that field-equivalent lift, but both the scalar
ConvertFrom path and the R2 M2-to-M1 expansion assume a bounded nonnegative
representative. A direct import/export and multiplication by one both return
the wrong field value.

CatBLST reproduces this in the focused portable-Debug test
`AvxMatrixNoKPrimitives.CanonicalImportKeepsBothRnsBasesConsistent`.

### Proposed repair

Generate the eight ConvertTo fixed-point quotient coefficients with floor
rather than ceiling. Leave R1, R2, and ConvertFrom unchanged.

For canonical radix-`2^50` inputs:

- each floor coefficient is exactly one less than the current ceiling value;
- the aggregate estimate cannot exceed the exact CRT quotient;
- the total fixed-point error is less than one quotient unit, so the estimate
  is either exact or one low;
- the corrected representative is consequently in `[0, 2p)`, within the
  existing nonnegative R2/ConvertFrom design envelope.

This is a generated-table-only production change: it adds no instruction,
branch, load, temporary, or runtime abstraction.

### CatBLST repair

Commit `5f0399e` implements the repair in the generator and generated table.
Only ConvertTo uses floor estimates; R1, R2, and ConvertFrom retain their
ceiling estimates. Each of the eight changed BLS12-381 constants is exactly
one less than its inherited VROOM value.

The Lean packet independently regenerates the floor table and proves:

- the estimated quotient never exceeds the exact CRT quotient;
- it is either exact or one low for a canonical radix-`2^50` input;
- the repaired representative lies in `[0, 2p)` and represents the imported
  canonical field value; and
- the exact prior witness now selects `p - t`, lies inside the existing R2
  input range, and preserves the canonical field value.

The exact portable-Debug import/export and multiply-by-one regression is now
enabled and passes. The constant validator and generated-Lean source check
pass, the complete Lean target builds without `sorry`, `admit`, axioms, or
`native_decide`, and 17 focused MatrixNoK tests pass. The remaining G1 FFT
test stops earlier in scalar `Fr::redc` Debug validation and is tracked
separately as `AUD-222`; it does not reach the repaired MatrixNoK comparison.

### Remaining evidence before upstream submission

- Resolve or bypass `AUD-222` with a valid scalar FFT oracle, then pass that
  remaining focused differential.
- Pass the repaired packet on native AVX-512IFMA.
- Confirm optimized import code differs only in embedded constants.
- Prepare the VROOM PR against its generator and tests, not only its checked-in
  generated BLS12-381 table.
