# Upstream Review Wave 01

## Record

- Local commit: `67a2dfcadf75f5452d503acc41ab3ebbb4568a5d`
- Upstream head: `54e6e55674722fc2797ebb4bbb71b26d881eb4b8`
- Review mode: parallel read-only comparison
- Benchmarks: not reviewed or run
- Code changes: none

This wave prioritized commits with likely safety, correctness, resource, ABI,
and performance consequences. A `port` or `translate` decision remains open
until implementation and validation are linked from the canonical manifest.

## C and C++ API lane

### `e61acc21931e4a8ef4ea924f97e6a44e6af154f2`

Decision: `not_applicable`.

At wave review, local `bindings/blst.hpp` still had the upstream defects. The
unsupported C-over-ABI wrapper, its commented test, SWIG C++ branch, Rust
package entry, and public documentation were removed in Phase 6A rather than
porting wrapper fixes. C++ support continues to refer to the native C++
implementation and API.

Required validation:

- retained C, native C++, and Rust builds do not include the wrapper;
- no package or install rule publishes it;
- no retained documentation advertises it; and
- repository search finds only historical audit references.

Residual risk: accidental continued packaging makes an unsupported, defective
API appear supported until deletion lands.

### `1519e46092a1fdef235eec79247c37f2c950fff0`

Decision: `port`.

`blst/blst_t.hpp:517-522` implements both `operator>>=` and `operator>>` using
`lshift_mod_256`. This is the exact defect fixed upstream. Replace both calls
with `rshift_mod_256` and add direct tests; no current test names this type or
operator.

Performance impact: none beyond selecting the correct existing primitive.

Residual risk: audit callers to determine whether prior results could have
escaped through a production API.

### `cf739f07477086e69f73b6f658c003249c356af4`

Decision: `translate`.

The arbitrary-length modular imports `blst_fp_from_le_bytes` and
`blst_fp_from_be_bytes` are absent. The fork uses split public headers and C++
exports, so a direct cherry-pick is inappropriate. Add the API to
`bindings/blst_scalar.h`, implement it in the local export architecture, and
define zero-length/null behavior deliberately.

Required validation:

- lengths 0, 1, 47, 48, 49, and multiple 48-byte chunks;
- both endian orders;
- differential reduction against an independent big-integer oracle;
- null behavior consistent with the C API contract.

Performance impact: linear in input chunks and outside core field kernels.

### `d4ca0f29b86c8512c786c8f50c5ffb5a9b1b6aa6`

Decision: `not_applicable`.

This applied only to `bindings/blst.hpp`. Phase 6A removed that unsupported
wrapper instead of maintaining its compiler-version policy.

### `f262a6e9985f84e1d2842960a158dc768b217884`

Decision: `translate`.

The relevant logic moved from `blst.h` to `bindings/blst_types.h`. It still
directly defines `bool` as `_Bool` for C99 through transitional C2x modes.
Translate upstream's use of `<stdbool.h>` for ordinary C, while limiting
special parser accommodations to retained Rust generation requirements. Go and
Zig compatibility must not be retained accidentally.

Required validation: standalone C99, C11, transitional C2x, final C23, C++,
and Rust-bindgen header parsing.

## Scalar and API hardening lane

### `dd5e80be97d03ef0e3568be4319f95434ceb8274`

Decision: `already_equivalent`; status `closed`.

The local implementation independently fixed the zero-length underflow in
commit `d6a55c08`. `blst/exports.cpp:410` and `:443` clear output and return
false before any pointer arithmetic for little- and big-endian imports. Direct
C API coverage at `test/C_api_tests.cpp:1809` includes `nullptr, 0` and output
clearing; header-level coverage is at `test/BindingHeader_tests.cpp:174`. Both
tests passed under `x64-Debug` during this review.

Performance impact: one predictable public-boundary branch.

Residual risk: low. The stronger output-clearing behavior should be documented,
and a Rust empty-slice regression remains desirable.

### `a3c15d6a6c6814d4d0d264a772c62fb86d4d41fe`

Decision: `port`, high priority.

`blst/sign/keygen.cpp:46-48` calls `blst_sha256_bcopy_v` for every salt length
up to 64, including zero. Portable copy happens to be safe; generated x86-64
and ARM64 assembly copies are do-while loops that copy once and decrement zero
to `SIZE_MAX`. Rust `key_gen_v5` passes the non-null pointer returned by an
empty slice and reaches this path. The existing C test checks `NULL, 0`
rejection, not a valid non-null empty salt.

Required validation:

- non-null zero-length salt in C and Rust;
- deterministic expected output;
- portable and assembly configurations;
- existing null-pointer policy remains unchanged.

Performance impact: a single initialization-time branch; no effect on nonempty
salt processing beyond that branch.

Residual risk: the assembly copy primitive retains a nonzero-length
precondition that must be documented and checked at all callers.

### `84926f21eb2bf15035a931ebab6481272290f251`

Decision: `translate`.

The fork independently replaced unbounded stack allocation in
`blst/scalar/curve/multi_scalar.cpp:272`, but allocates the complete remaining
table on the heap at `:323-337`. Allocation failure silently zeroes a caller
buffer through a void API. Preserve the fork's checked size arithmetic while
adapting upstream's bounded, allocation-free chunking.

Owner decision: the production implementation must be allocation-free; do not
retain a heap fallback for oversized precomputation.

Required validation:

- G1 and G2 across window sizes, including threshold crossings;
- packed and pointer-array inputs;
- infinity points and table reuse;
- an intentionally low scratch limit to force chunking;
- evidence that the production path performs no heap allocation.

Performance impact: removes nondeterministic allocation latency. The exceptional
large case may perform more conversions and must be measured later with only
one benchmark process at a time.

## Rust safety lane

### `1172c188cfcc059121b024a2010884700402baf9`

Decision: `translate`, high priority.

The native pairing context retains a DST pointer. Local
`bindings/rust/src/pairing.rs` exposes a non-generic safe `Pairing`, accepts an
unconstrained `&[u8]` in `new` and `init`, and stores no lifetime owner. A
caller can drop temporary DST storage and later cause native code to read freed
memory.

Implementation decision:

- accepted design: own DST in trailing storage within the existing aligned
  allocation. This preserves the non-generic Rust API and uses one allocation.
  Reinitialization reuses sufficient capacity or allocates and swaps safely.

Required validation: compile-fail lifetime tests for the borrowed design, or
temporary-vector lifetime/Miri tests for the owned design; reinitialization,
empty DST, threading, and merge behavior.

Residual risk: high until fixed. Audit the local unsafe `Send` implementation
and native ownership documentation with the selected design.

### `96a46e6ca7360713d9bb1fb78469c28b3c043099`

Decision: `translate`, high priority.

Local `Pairing::mul_n_aggregate` passes `scalar.as_ptr()` and caller-provided
`nbits` directly to native G1/G2 aggregation without proving the slice holds
the requested bits. Use an overflow-safe ceiling:

```text
nbits / 8 + (nbits % 8 != 0)
```

and reject undersized slices before type dispatch or FFI.

Required validation: undersized, exact byte-boundary, zero, extreme `nbits`,
both G1/G2 branches, and high-level randomized aggregate APIs.

Performance impact: one comparison outside expensive pairing work.

Residual risk: audit other pointer/bit-length pairs and planner size arithmetic.

### `62719185360d9e0887bb5a1ba68c6b89d61b01c9`

Decision: `translate`, high priority.

The local threaded MSM scheduler uses its small path only below 32 points, then
calls `breakdown` even when worker count greatly exceeds point count.
`breakdown` can reduce its window divisor to zero and divide by it. Translate
upstream's oversubscription routing and all planner guards. Because the fork
still uses an unbounded result channel, assess upstream's low-overhead
Barrier/Cell accumulation as a coherent performance change rather than
mechanically applying fragments.

Owner decision: implement the routing and planner guards, and retain the
upstream low-overhead per-worker result-slot direction subject to its unsafe
initialization and barrier/panic audit.

Required validation: direct planner properties over bit/window/CPU ranges,
injectable oversubscribed worker counts, G1/G2 differential MSM, and relevant
feature combinations.

Residual risk: medium-high availability risk until fixed; after translation,
review uninitialized storage and worker-panic/barrier behavior.

## Assembly and generated-output lane

### Trust-chain prerequisite

Decision: repair before porting assembly commits.

`build/refresh.sh` still searches `../src/asm`, although generators now live
under `blst/asm`. It also couples assembly generation to definition-file and
Rust-bindgen regeneration. Checked-in Win64 outputs contain manual
`BLST_PORTABLE` changes while their generators still emit
`__BLST_PORTABLE__`. The repository therefore cannot currently attest that
checked-in assembly came from the checked-in generators.

Required validation:

- split or narrowly scope the assembly generator;
- regenerate twice with byte-identical second output;
- enforce an expected output allowlist;
- reconcile generator and MASM feature-macro policy.

### `525a87a` and `e789f8a`

Decision: `translate` both as generator change plus generated companion.

Remove obsolete `sqr_mont_382x`. No production caller was found; local removal
also includes fork-specific declarations/wrappers and its interface test.
Unsupported internal symbol consumers are the principal compatibility risk.

### `4a74ffd` and `e789f8a`

Decision: `port` the generator schedule; `translate` the generated refresh.

Local x86-64 generators retain the older ADX carry chains for sparse-256,
Montgomery-384, and non-reduced-383 routines. Preserve the upstream schedule
exactly. It affects scalar Fr/Fp directly and AVX-512 indirectly wherever AVX
APIs delegate serialization, validation, signing, or conversion to scalar.

Required validation: ADX/non-ADX/portable differential arithmetic, modulus
boundaries, aliasing permutations, and scalar/AVX uncompression paths. Object
code must match the intended instruction schedule before performance evidence
is considered.

### `e3d708b`

Decision: `port`.

Update the x86 translator so `{nf}` prefixes survive NASM translation. No
current generator uses the prefix, making this tooling/future compatibility
rather than an immediate production-code change.

### `33d7f80` and `bf008ec`

Decision: `port` the corrected contract and `translate` the call-site changes.

The tailored sparse-256 multiply permits a non-reduced operand only in the
second multiplicand position. Local `Fr::to_mont` currently calls multiplication
with possibly unreduced input first and the reduced constant second. Although
ordinary multiplication is commutative, the implementation's operand bounds
are not symmetric. Swap the operands and audit every related conversion,
reciprocal, parsing, secret-key, and AVX FFT helper path.

Required validation: portable/assembly differential tests across reduced and
unreduced boundaries, both aliasing directions, and conversion round trips.

### `41b511a` and `14cd09a`

Decision: `translate` both as generator policy plus generated companion.

Local generators and Win64 outputs lack guards preventing fallback `$1`
symbols from being exported or referenced in non-portable MSVC objects. The
CMake Win64 build includes all x86 assembly objects. Define how the portable
feature macro reaches MASM, then apply the guards and regenerate.

Required validation: MSVC/MASM portable and non-portable builds plus
`dumpbin /symbols` or equivalent proving fallback symbols are absent or
resolved as intended.

## Consolidated implementation order

1. Fix `AUD-004` empty-salt assembly-path crash.
2. Fix `AUD-001` and `AUD-002` Rust safe-API memory hazards.
3. Fix `AUD-008` Fr operand-order correctness gap.
4. Fix `AUD-003` threaded MSM divide-by-zero.
5. Repair deterministic assembly generation (`AUD-009`) and MASM macro policy.
6. Port the accepted generator chain and regenerate once.
7. Translate bounded MSM precomputation (`AUD-005`).
8. Port the retained C/C++ correctness and API changes.
