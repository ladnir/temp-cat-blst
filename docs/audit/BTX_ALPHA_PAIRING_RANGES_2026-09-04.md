# BTX Alpha Pairing Range Review

## Purpose

After BTX verifies a combined G1 decryption share, it computes one target-group
alpha value for every admitted ciphertext. Every alpha uses the same combined
G1 point and one contiguous entry from the epoch's immutable prepared-G2
table. The existing prefix interface is allocation-free, but a parallel caller
would need one full-batch AVX workspace per worker. That unnecessarily
multiplies resident scratch by the batch size.

This change adds contiguous scalar and native-AVX pairing ranges. The AVX
workspace constructor accepts the maximum range size, so one worker retains
only the Miller and final-exponentiation scratch required by its grain. The
interfaces are generic prepared-pairing operations; CatBLST contains no BTX
task, priority, cancellation, or phase state.

## Interface and ownership

The scalar C operation is
`blst_scalar_prepared_g2_batch_pairing_range`. It accepts an immutable prepared
batch, one library-produced scalar G1 projective point, a checked `[begin,
begin + count)` range, and caller-owned GT output. It converts the common G1
point to affine once, evaluates each prepared line table, and final
exponentiates each result. It allocates no memory and needs no mutable
workspace.

The AVX additions are
`blst_avx_prepared_g2_batch_workspace_new_with_capacity` and
`blst_avx_prepared_g2_batch_pairing_range`. A workspace owns aligned Miller
slots and batch-final-exponentiation scratch for at most `capacity` outputs.
One range fills those slots from consecutive prepared operands and invokes the
existing batch final-exponentiation kernel. The prepared table may be shared
across workers. A mutable workspace and output range may not be shared.

Both Rust namespaces expose initialized-output methods for ordinary callers
and `pairing_range_uninit_into` for task managers. The latter accepts
`MaybeUninit<Gt>` storage. Success establishes every output's native GT
invariant; after failure, the caller must continue treating the complete range
as uninitialized. This lets BTX allocate output capacity without first writing
identities that the native kernel immediately overwrites.

## Preconditions and failure behavior

- Prepared G2 and G1 values are library-produced, subgroup-validated values.
- The contiguous range must fit inside the immutable prepared batch.
- AVX output count must not exceed workspace capacity.
- Nonempty output must not overlap the G1 input, prepared batch, or AVX
  workspace.
- Concurrent calls require disjoint outputs and, for AVX, distinct workspaces.
- Timing may depend on the public range length and public point identities.
- Empty ranges succeed and publish nothing.
- Invalid pointers, ranges, or capacities are rejected before output writes.
- No C++ exception crosses the C boundary. If an unexpected exception occurs
  after arithmetic begins, output is unspecified; the Rust uninitialized
  method therefore establishes initialization only on success.

The capacity constructor rejects zero and capacities larger than the prepared
batch. BTX creates no workspace for an empty alpha phase. This keeps the
low-level contract unambiguous without requiring a zero-capacity native owner.

## Memory and performance properties

Range execution performs no allocation. Scalar stack usage is constant. An
AVX workspace retains `capacity` Miller values plus fixed batch-final-
exponentiation scratch, rather than storage for the complete epoch table.
Prepared lines remain immutable and are not copied into the workspace.

The Monad task currently uses a provisional grain of 16. That value is not a
settled performance policy. A later sequential Peach A/B must compare the
existing full-prefix call with serial range calls and then measure 2-, 4-, and
pool-width execution. No benchmark was run in this correctness tranche, and
no two benchmark processes may run concurrently.

## Correctness evidence

The scalar Rust and C tests compare split contiguous ranges with direct
pairings and reject out-of-range requests without changing canaries. The AVX
Rust and C tests compare a bounded-workspace subrange with direct pairings and
cover zero, oversized, and out-of-range capacities and ranges. The existing
initialized-output methods delegate to the same checked uninitialized-output
boundary.

On Peach, Clang 21 release/native builds pass the focused scalar and AVX C API
tests. Exact-source CatBLST Rust runs pass 75 scalar and 76 AVX-feature library
tests. Monad exercises both uninitialized-output wrappers through its complete
scalar and native-AVX suites: 120 scalar unit tests, 143 AVX-enabled unit
tests, and seven conformance tests in each configuration pass. Dedicated
Clang 21 ASan/UBSan builds also pass the scalar-default C API test and the
native-AVX compatibility C API test with leak detection and halt-on-error
enabled. The scalar test is intentionally run in the scalar-default executable
because the public default GT representation follows the configured backend.

## Review conclusion

The interface preserves CatBLST's kernel boundary and BTX's native-type
invariant. It avoids full-batch scratch replication and avoids redundant GT
initialization. Measured grain selection remains open; it does not require a
protocol or caller-facing API change.
