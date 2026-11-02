# Middle Product

The Rust middle-product surface is intentionally small:

- `MiddleProductPrecomputedB`
- `precompute_middle_product_b`
- `middle_product`
- `update`
- `combine`
- `ScalarMiddleProductPrecomputedB` with `middle_product_scalar` and
  `update_scalar`, plus `ScalarMiddleProductWorkspace` and
  `middle_product_scalar_with_workspace`
- `AvxMiddleProductPrecomputedB` with `middle_product_avx` and `update_avx`
  plus `AvxMiddleProductWorkspace` and
  `middle_product_avx_with_workspace`
- `NativeAvxMiddleProductPrecomputedB` with
  `NativeAvxMiddleProductWorkspace` and
  `middle_product_avx_native_with_workspace`

Design points:

- the compatibility and explicit scalar interfaces take `b` as scalar-ABI
  `&[blst_p2]`, `a` as `&[blst_p1]`, and `a0` as `&blst_p1`
- those interfaces write output into caller-owned `&mut [blst_fp12]`
- precomputation owns an opaque native handle
- compatibility evaluation and update calls allocate native staging/work
  buffers on each call
- the explicit `*_with_workspace` evaluations reuse caller-owned native FFT,
  output-staging, and AVX input-conversion storage and allocate nothing after
  successful workspace construction
- the native AVX interface accepts `blst_avx_p2` and `blst_avx_p1` directly,
  writes caller-owned `blst_avx_gt`, and needs neither conversion storage nor
  scalar output staging

## Stateless range execution

The ranged interface exposes the internal middle-product schedule without
putting a scheduler in CatBLST. It is intended for a caller that already owns a
thread pool, task priorities, cancellation policy, and phase barriers.

The scalar Rust types are `ScalarMiddleProductPrecomputedB` and
`ScalarMiddleProductWorkspace`. Their range methods are:

- `ScalarMiddleProductPrecomputedB::range_work_count`;
- `ScalarMiddleProductWorkspace::input_range`;
- `ScalarMiddleProductWorkspace::execute_range`; and
- `ScalarMiddleProductWorkspace::output_range`.

The conversion-free AVX equivalents are
`NativeAvxMiddleProductPrecomputedB` and
`NativeAvxMiddleProductWorkspace`. They expose the same four method names and
use native `avx::G1Projective` and `avx::Gt` values throughout.

The corresponding C families are
`blst_middle_product_scalar_range_work_count`,
`blst_middle_product_scalar_input_range`,
`blst_middle_product_scalar_execute_range`, and
`blst_middle_product_scalar_output_range`, with
`blst_middle_product_avx_native_*` counterparts. They are declared in
`bindings/blst_middle_product.h`. The Rust methods retain the same work domains
and safety contract.

### Required execution order

One evaluation has the following state machine:

| Order | `MiddleProductRangePhase` | Stage | Barrier |
| ---: | --- | ---: | --- |
| 1 | `InputLoad` | 0 | After every input range. |
| 2 | `ForwardBitReverse` | 0 | After every swap range. |
| 3a | `ForwardJ0` | each forward FFT stage | Before `ForwardNonzero` for that stage. |
| 3b | `ForwardNonzero` | the same stage | Before the next forward stage. |
| 4 | `AffineConversion` | 0 | Before pointwise pairing. |
| 5 | `PointwisePairing` | 0 | Before the inverse transform. |
| 6 | `InverseBitReverse` | 0 | Before inverse butterflies. |
| 7 | `InverseButterflies` | each inverse FFT stage | Before the next inverse stage. |
| 8 | `Normalization` | 0 | Before extraction. |
| 9 | `Extraction` | 0 | Before the caller reads any output. |

For one row, query its work count, partition `0..count` into disjoint ranges,
run every range, and wait for all of them. Range completion order within that
row is irrelevant. The caller must not overlap two different phases or stages.

`InputLoad` and `Extraction` use their typed methods. Every phase between them
uses `execute_range`. The input domain covers the complete padded FFT buffer:
indices below the logical input length copy a point and the remaining indices
write the identity. The extraction domain covers the logical output length.

Forward and inverse stage counts equal `log2(fft_work_count)`, where
`fft_work_count` is the `InputLoad` count. A caller may instead cache all valid
stage counts when it constructs its immutable plan. It should not discover
stages repeatedly in a hot scheduling loop.

### Ownership and concurrency

The precompute is immutable after construction and may be shared by every
range. The workspace owns the mutable, 64-byte-aligned FFT buffers. Its address
and capacity remain stable for its lifetime. A workspace constructed for a
larger FFT capacity may execute a smaller matching precompute by using the
required prefix.

Concurrent calls may share one workspace only when all calls address disjoint
ranges in the same phase and stage. This deliberate shared-mutable boundary is
`unsafe` in Rust. The caller must additionally keep the precompute, input,
output, and workspace alive until all submitted ranges complete. Dropping a
future, cancelling a logical task, or returning from the scheduling function
does not end the lifetime of a range that is already running.

Each affine-conversion range performs one inversion. Excessively fine affine
partitioning can therefore lose performance even when other phases benefit
from the same grain. Choose range grain per phase from end-to-end
measurements.

### Failure and publication

Empty valid ranges are successful no-ops. An invalid handle, phase, stage,
range, length, backend, or insufficient workspace capacity returns an error
before the addressed operation writes. `output_range` treats its destination
as uninitialized writable storage; it initializes every addressed element on
success and does not read its previous value.

CatBLST does not retain a phase cursor. It therefore cannot reject a
well-shaped range merely because the caller scheduled it in the wrong order.
Missing a barrier, concurrently overlapping ranges, or executing phases out of
order violates the API contract and leaves the workspace result unspecified.
The caller must publish output only after every extraction range succeeds.

Range calls allocate nothing and create no threads. They are variable-time
with respect to public points, sizes, ranges, and schedules. The one-shot
`middle_product_scalar_with_workspace` and
`middle_product_avx_native_with_workspace` operations remain useful sequential
correctness oracles for a new ranged driver.

### Minimal scheduler shape

The following pseudocode shows the required control structure. `run_phase`
partitions one work domain, submits disjoint ranges, and returns only after all
ranges succeed.

```rust,ignore
workspace.input_ranges_then_barrier(input, plan)?;
run_phase(ForwardBitReverse, 0)?;

for stage in 0..plan.forward_stage_count() {
    run_phase(ForwardJ0, stage)?;
    run_phase(ForwardNonzero, stage)?;
}

run_phase(AffineConversion, 0)?;
run_phase(PointwisePairing, 0)?;
run_phase(InverseBitReverse, 0)?;

for stage in 0..plan.inverse_stage_count() {
    run_phase(InverseButterflies, stage)?;
}

run_phase(Normalization, 0)?;
workspace.output_ranges_then_barrier(plan, output)?;
```

The application may implement `run_phase` with Rayon, a dedicated worker
pool, or direct serial calls. Calling one full range per phase preserves the
same state machine and provides a near-zero-overhead sequential migration
path.

Mathematical and point-domain contract:

- for `a.len() == n` and `b.len() == 2*n`, output `k` is
  `product_i pairing(a[i], b[n + k - i])`
- every `blst_p1` and `blst_p2` input must already be reduced, on curve, and in
  its prime-order subgroup; infinity is accepted as the identity
- the native implementation performs no point or subgroup validation
- outputs are GT elements stored in `blst_fp12`; `combine` likewise requires GT
  inputs, not arbitrary Fp12 values

Length rules:

- `precompute_middle_product_b` expects `b.len() == 2 * n`
- `middle_product` expects `a.len() == n`
- `update` writes exactly `n` outputs for the precomputed handle
- `combine` expects `mid.len() == cnst.len()`
- `out.len()` must be at least `n`
- unsupported native FFT lengths fail precomputation rather than wrapping size
  arithmetic or attempting a transform

Backend behavior:

- compatibility interfaces keep the scalar/public Rust ABI
- the default handle follows the CatBLST build default for compatibility
- the explicit scalar and AVX handles bind the backend at construction and
  have different Rust and C types; they cannot be interchanged
- the native AVX handle and workspace never cross the scalar representation;
  canonical encoding and decoding remain caller-owned wire-boundary work
- an explicit AVX precompute returns `PrecomputeFailed` when the linked build
  does not contain AVX support
- the AVX implementation converts from scalar boundary types once per kernel
  call into its caller-owned conversion buffer; no per-point backend tag,
  allocation, or inner-loop dispatch is present
- `combine` uses the GT group law; at the `blst_fp12` boundary this appears as
  elementwise multiplication because GT is represented inside Fp12
- all algorithms are variable-time with respect to their public points and
  lengths
- the opaque handle is immutable after successful construction; the safe Rust
  wrapper owns its lifetime and does not expose concurrent mutation
- workspaces are mutable and not shared by concurrent evaluations; callers
  needing independent complete evaluations create one workspace per
  evaluation; ranged execution may deliberately share one workspace under the
  stricter disjoint same-phase-and-stage contract above
