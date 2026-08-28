
# Lean typed AVX curve API

The safe `catblst::avx` types retain G1, G2, and GT values in CatBLST's native
AVX representation. They are intentionally distinct from the corresponding
`catblst::scalar` types: the representations have different sizes, alignment,
and internal layouts. A protocol chooses a backend before constructing these
values and does not convert resident arithmetic values between backends.

The module also retains its raw `blst_avx_*` expert ABI. New Rust protocol code
should use the typed layer unless it is implementing another reviewed
allocation-free primitive inside CatBLST.

## Invariants

Safe typed values carry these invariants:

- `G1Projective` is an on-curve point in the prime-order G1 subgroup;
- `G2Projective` is an on-curve point in the prime-order G2 subgroup;
- `Gt` is an element of the prime-order pairing target group; and
- the identity is valid for every group type.

Canonical decoders establish the curve and subgroup invariants before
returning. Safe arithmetic preserves them and therefore does not repeat
validation in each kernel. Protocol admission, authorization, and proof
verification remain caller responsibilities.

The native byte arrays inside raw AVX structs are opaque arithmetic storage.
They are not serializations and must not be inspected, modified, hashed, or
sent over a wire. Safe typed wrappers do not expose mutable raw storage.

## Supported typed surface

The lean surface provides only the operations required by current native
consumers:

- G1 and G2 identity, generator, addition, fixed-width scalar multiplication,
  canonical compressed encoding, strict decoding, equality, and identity
  testing;
- GT identity, multiplication, squaring, inversion, fixed-width scalar
  multiplication, pairing, canonical encoding, strict decoding, and equality;
- immutable prepared-G2 line tables;
- a fixed-capacity batched-pairing workspace;
- caller-owned scratch for randomized G2 batch subgroup validation;
- a caller-owned projective G1 MSM workspace; and
- typed native middle-product precomputation and evaluation.

This is not a source-compatible AVX clone of `blstrs` or of the scalar module.
Additional convenience operations should be added only when a real consumer
requires them and their allocation, timing, validation, and aliasing contracts
are understood.

## Encoding and backend compatibility

G1 and G2 use CatBLST's canonical compressed encodings of 48 and 96 bytes.
GT uses twelve canonical 48-byte big-endian Fp coefficients, ordered by Fp2
index, then Fp6 half, then Fp component. The scalar and AVX encoders produce
identical bytes, and either backend may decode bytes produced by the other.

Canonical encoding is the supported boundary for serialization, hashing,
cross-backend conformance tests, and backend selection. It is not an internal
arithmetic strategy. Typed AVX operations never convert through scalar point
or Fp12 storage.

## Timing

Safe G1, G2, and GT scalar multiplication accepts `crate::Scalar`, so the
canonical-scalar precondition is represented by the type. Each operation uses
a fixed public 255-bit width; runtime and memory access are independent of the
scalar contents. These operations allocate nothing and securely erase their
temporary canonical byte encoding. GT scalar multiplication calls CatBLST's
constant-time cyclotomic exponentiation kernel directly; it does not branch in
Rust on scalar bits.

Decoding, subgroup validation, encoding, equality, identity tests, affine
conversion, pairing, final exponentiation, and variable-time MSM kernels do
not carry a constant-time promise. Their inputs must be public unless a more
specific native contract states otherwise.

`G1MsmWorkspace::msm_into` accepts native typed `Scalar` values and requires
them to be public. It uses shared-doubling GLV wNAF-5 batches through the
measured small/medium crossover and Pippenger above it. Both algorithms branch
and address tables using scalar digits. BTX uses this method for public
Lagrange coefficients during share combination. The low-level
`pippenger_msm_into` method retains the canonical-byte Pippenger path for
explicit comparison and specialized callers.

`G1MsmWorkspace::secret_msm_into` accepts native typed `Scalar` values and uses
a constant-time GLV MSM for secret coefficients. Its control flow and memory
addresses depend on the public point count, but not on scalar contents. The
kernel shares its doubling schedule in 32-, 16-, and 8-point chunks, then uses
fixed-width constant-time multiplication for the short tail.

## Allocation and workspaces

Single point, GT, and pairing operations use stack or caller-provided output
storage and allocate nothing. `PreparedG2Lines::new`, `PreparedG2Batch::new`,
`PreparedG2BatchWorkspace::new`, `PairingWorkspace::for_pairs`,
`G1MsmWorkspace::for_points`, and native middle-product constructors allocate
reusable storage outside the hot path. They report allocation failure.

`avx_p2s_in_g2` requires a mutable byte slice of at least
`blst_avx_p2s_in_g2_workspace_sizeof()` bytes. The native validation kernel
does not allocate. The safe wrapper may allocate its temporary pointer list;
the infrequent public-parameter admission path is not a resident arithmetic
hot loop.

After construction:

- `PreparedG2Lines::pairing` allocates nothing;
- `PreparedG2BatchWorkspace::pairing_prefix_into` and
  `pairing_indices_into` allocate nothing;
- `PairingWorkspace::clear`, `push`, and `evaluate` allocate nothing;
- `G1MsmWorkspace::msm_into` and `secret_msm_into` allocate nothing; and
- typed native middle-product evaluation allocates nothing.

`PairingWorkspace::for_pairs` verifies backend availability and records the
requested logical pair limit rather than relying on allocator capacity.
`PairingWorkspace::push` converts both operands before mutating either staging
array and fails before mutation when the logical capacity is exhausted.
All G1 MSM methods check lengths and scratch capacity before calling the native
kernel. Variable-time Pippenger initializes every bucket before its first read.
The public wNAF path completely overwrites every live digit and odd-multiple
table entry; those values are public and are not pointlessly erased. Secret MSM
stores its public point tables and split-scalar staging in the caller-owned
workspace rather than on the thread stack. It erases the split scalars,
canonical scalar bytes, and Booth indices before returning; public point tables
may remain. Callers use one mutable workspace per concurrent operation.
Immutable prepared line tables and native middle-product precomputations may be
shared after construction.

`PreparedG2Batch` owns one 68-line native Miller table and one identity flag
per G2 operand. Construction default-initializes the large table allocation
without clearing it, then completely overwrites every nonidentity table before
publication. Identity entries are never read as line tables. The
`storage_bytes` result accounts for the object, tables, and flags but not
allocator metadata. At the current AVX ABI, line storage is approximately
68 KiB per operand; callers must choose the prepared batch extent as a public
parameter and account for that resident memory.

`PreparedG2BatchWorkspace` owns one Miller-result slot per prepared operand and
establishes the capacity and logical size of all six structure-of-arrays
buffers used by batch final exponentiation. Thus even the first evaluation
does not grow or initialize those buffers. Prefix evaluation maps output
position `i` to prepared operand
`i`. Indexed evaluation maps output position `i` to the public index supplied
at `indices[i]`; repeated indices are permitted. Both operations convert the
shared G1 operand to affine once, test its identity status once, run only the
finite-point Miller inner loop for nonidentity pairs, and batch final
exponentiation over the requested result count. A workspace may be reused with
another prepared batch only when its recorded capacity covers the request.
One immutable batch may be shared, but each concurrent evaluation requires a
distinct mutable workspace and output slice.

## Failure and validation

Safe methods return `catblst::avx::Error`; they do not panic for malformed
external input or backend unavailability. Strict point decoding distinguishes
malformed encodings from subgroup failure. GT decoding rejects noncanonical Fp
coefficients and values outside GT while accepting the identity.

The safe layer checks shape, capacity, and backend availability where those
conditions can arise from ordinary caller input. Native kernels may still rely
on the safe-type invariants and fixed aliasing established by Rust borrows.
Debug native builds assert additional algebraic preconditions.

Raw prepared-G2 construction reports backend unavailability, invalid native
layout, allocation failure, and initialization failure separately. The safe
prepared-line wrapper preserves unavailability and allocation errors and maps
native layout or initialization rejection to `Error::BackendFailure`.
Prepared-batch construction reports `EmptyPreparedBatch` for empty input and
reports native allocation or preparation failure separately. The safe
evaluators report pairing length, workspace capacity, and out-of-range index
errors before entering the native kernel.
The C boundary likewise validates all sparse indices before publishing output;
after a backend exception, output contents are unspecified.

## Deliberate omissions

The typed AVX layer does not provide dynamic backend dispatch, per-value
backend tags, implicit scalar/AVX conversions, hidden hot-loop allocation,
thread pools, Rayon integration, `ff`/`group`/`pairing` traits, Serde support,
or a complete mirror of the scalar convenience API. Higher-level code may use
type erasure only outside arithmetic kernels after fixing all associated
backend types.
