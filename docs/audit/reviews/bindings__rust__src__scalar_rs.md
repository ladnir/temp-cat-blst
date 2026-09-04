# File Review: `bindings/rust/src/scalar.rs`

## Record

| Field | Value |
| --- | --- |
| Path | `bindings/rust/src/scalar.rs` |
| Review scope | lean native scalar curve surface and prepared-G2 range addition |
| Reviewed base | `295c17c4` plus the uncommitted 2026-09-03 candidate |
| Date | 2026-09-03 |
| Disposition | in-scope production candidate, pending owner approval |
| Profiles | scalar Rust/C ABI; selected AVX types remain separate |
| Related files | `bindings/blst_scalar.h`, `bindings/rust/src/bindings.rs`, `bindings/rust/docs/scalar.md`, `blst/c_api.cpp` |

## Purpose and trust boundary

This module provides the lean Rust scalar-field, G1, G2, Miller-result, and GT
types used by protocol code that selects the scalar backend. Safe constructors
establish reduced-field or prime-subgroup invariants. Safe arithmetic preserves
those invariants and intentionally does not repeat curve or subgroup checks.
Untrusted bytes cross only explicit decoding functions. Raw mutation and raw
construction remain unsafe because they can violate every later arithmetic
precondition.

The prepared-G2 addition moves repeated public G2 line preparation outside a
hot pairing-product operation. It does not widen the accepted algebraic domain
or add an alternate protocol flow.

## Representation and ownership

Scalar, point, affine-point, Miller, and GT wrappers are transparent over the
matching scalar C storage. They do not allocate or own external resources.
`PreparedG2Batch` is the explicit exception: it owns one opaque native heap
object through `NonNull`, caches its validated length and storage byte count,
and frees it exactly once through `Drop`.

The owner is neither cloneable nor serializable. It is `Send + Sync` because
the C object is completely initialized before publication, evaluation accepts
only a const pointer, and all per-call state is on the calling thread's stack.
Rust ownership prevents destruction from racing with a safe borrow.

## Prepared-G2 contracts

- `new` rejects an empty slice before FFI. Its safe G2 inputs already establish
  on-curve and subgroup membership; identity is valid.
- Native construction returns either a complete owner or null. Rust maps null
  to `AllocationFailed`, then cross-checks the native length and nonzero storage
  report before publishing the safe value.
- `miller_loop_product_indices` requires equally sized G1 and index slices. It
  locates and reports the first out-of-range index before FFI.
- Successful empty evaluation returns the pre-final multiplicative identity.
- Successful nonempty evaluation returns `MillerLoopResult`, not `Gt`. The only
  route from that type to `Gt` is explicit final exponentiation.
- The operation allocates nothing. Timing may depend on public length, index
  values, and identity locations.
- Independent calls may share the owner. Range products may be multiplied in
  any tree shape before one final exponentiation.

## Failure, memory, and concurrency analysis

Rust validates all slice relationships that it can express. The raw C layer
revalidates indices before arithmetic, so a wrapper bug cannot publish a
partial result. `MaybeUninit` becomes initialized only after C reports success.
Constructor postcondition checks free the native object before returning
`BackendFailure`; no partially checked owner escapes.

The unsafe surface is narrow: the native constructor receives a contiguous
transparent G2 slice, evaluation receives a contiguous transparent G1 slice
and initialized `usize` indices, and `Drop` consumes the unique native owner.
There is no aliasing mutable reference, shared workspace, callback, or lifetime
borrow hidden behind the safe API.

Construction and identity filtering are variable-time over public data.
Prepared values and points are not secret. The range operation uses bounded
stack arrays and no dynamic dispatch. It preserves scalar-native storage and
does not serialize or convert through AVX.

## API-surface decision

The module exposes one prepared owner and one indexed product operation. It
does not add separate prefix, range, pairing, or final-GT variants. Ordinary
slices express a range. Indexing supports reuse of BTX's reverse-ordered
prepared table without a second owner. Callers explicitly combine
`MillerLoopResult` values and final-exponentiate once.

## Validation

- Rust release/native: 73 library tests pass on Peach with Clang 21.
- Raw C++ release/native: the focused test passes, and the broad suite
  passes 580 tests with four unconfigured external-vector skips.
- ASan/UBSan: the focused prepared-batch test passes.
- Allocation failure: evaluation succeeds after all global C++ allocation
  routes are disabled.
- Differential tests cover direct pairing, arbitrary index permutations,
  identities, empty input, seven partition grains, reverse completion order,
  and concurrent immutable sharing.
- `check_rust_bindings.ps1` reports exact 258/258 C/raw-Rust declaration parity.
- `check_public_api_documentation.ps1` accounts for the family as `DOC-027`.

Detailed benchmark and command evidence is in
`docs/audit/BTX_SCALAR_PREPARED_PAIRING_RANGES_2026-09-03.md`.

## Residual risk and decision

The operation is suitable for the CatBLST scalar backend pending owner review.
The local range kernel does not select a Rayon grain or enforce a BTX phase
barrier; those are intentionally application responsibilities. The BTX
integration must use a parallel tree reduction and choose grain from logical
batch size and pool width. Cross-platform release matrices remain library-wide
evidence and are not required to establish this narrow Linux/Clang deployment
path.
