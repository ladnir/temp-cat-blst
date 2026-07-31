# Upstream Review Wave 02

## Record

- Local commit: `67a2dfcadf75f5452d503acc41ab3ebbb4568a5d`
- Upstream head: `54e6e55674722fc2797ebb4bbb71b26d881eb4b8`
- Review mode: parallel read-only comparison
- Commits reviewed: remaining 44
- Benchmarks: not reviewed or run
- Code changes: none

Wave 02 completed the decision pass for every upstream-only commit. Literal
Go, Zig, Java, Python, Node, Emscripten, and C# changes were closed only after
checking whether their defect class or test intent applied to retained
C/C++/Rust code.

## Rust retained surface

| Commit | Decision | Evidence and follow-up |
| --- | --- | --- |
| `07b109d` | translate | Local threaded Rust still uses unbounded result channels in aggregate verification, Miller loop, and Pippenger. Construct bounded `sync_channel` instances only after final worker counts; compose with `6271918`. |
| `b83544b` | port | `AggregatePublicKey::sub_aggregate` is absent; p1/p2 conditional-negation FFI already exists. Add the upstream additive API and G1/G2 group-identity tests. |
| `e7f90de` | reject | Upstream version `0.3.16` would falsely claim upstream release identity. CatBLST needs an independent crate name/version/publish policy. |
| `dece82e` | translate | Local dev manifest still directly depends on `rand 0.8`. Use `rand_core` traits with the retained ChaCha RNG, preserve MSRV, and compile-check all Rust targets without auditing benchmarks. |
| `23e330f` | not_applicable | Zig is removed. Its owned trailing-DST design strengthens the selected translation options for Rust `1172c18`; C++ already owns/copies DST. |
| `23cdd8a` | not_applicable | Go was removed in Phase 6B. The shared scalar-length bug remains actionable through Rust `96a46e6`; raw C/C++ pointer/bit-length contracts need documentation. |
| `54e6e55` | reject | Do not import upstream version `0.3.17`; resolve once under the CatBLST package identity. |

### Channel translation constraints

Upstream's bounded-channel optimization is not a mechanical substitution. Each
capacity must be derived after the actual worker count. The small Pippenger
channel disappears when `6271918` is translated, so it must not be
reintroduced. Add deadlock/time-bound stress tests around worker and stride
thresholds.

### Rust package identity

Publishing should remain disabled until the rename chooses:

- crate name and semantic-version starting point;
- repository and documentation links;
- MSRV and feature policy;
- lockfile/audit policy; and
- package dry-run and consumer migration checks.

## Retained C/C++/ARM/generated/support surface

| Commit | Decision | Evidence and follow-up |
| --- | --- | --- |
| `9c35e86` | translate | Generated Rust FFI companion for `cf739f0`; regenerate after native API implementation and tooling repair, never hand-edit or apply standalone. |
| `6aa1fe7` | not_applicable | Java/SWIG-only guarded byte wrapper and typemap changes; removed in Phase 6C. |
| `d7ab1f0` | not_applicable | The unsupported C-over-ABI wrapper and SWIG were removed in Phases 6A and 6C; native C++ API design remains separate. |
| `67fc54b` | defer | Reject Zig listing. Defer V6 conformance claim until retained protocol/vector evidence supports it. |
| `2e1948b` | port | Tracked CHERI generator/artifacts require `scbnds` translation support. |
| `ed865dc` | port | Bound 256/384-bit inverse scratch capabilities; depends on translator port and deterministic generation. |
| `fe85dc3` | translate | Mixed generated companion for CHERI inverse and Win64 `cf739f0` exports; regenerate from accepted source changes only. |
| `68b58bd` | translate | Move CodeQL to the current major, remove C#, add retained Rust, remove private-repo skip, and establish meaningful scalar/AVX C++ coverage. |
| `ff66582` | translate | Retain independent draft-v13 generator/compression/pairing vectors in C++ and Rust scalar/AVX tests; drop Python harness. |

### CHERI dependency chain

The checked-in CHERI sources and generated objects make the two source commits
applicable. If a later architecture-policy decision removes CHERI, remove its
generators and artifacts explicitly rather than leaving stale, unbounded
outputs.

Validation requires translator snapshots, CHERI assembly acceptance, inverse
vectors, capability-bound checks, ordinary ARM64 object-code equivalence, and
portable differential inversion tests.

### CodeQL translation

The current workflow is manual, skips private repositories, and analyzes C++
plus C#. C# is unsupported and Rust is retained. CodeQL is useful for C/C++ and
Rust source analysis but cannot validate handwritten/generated assembly,
constant-time behavior, or the scalar/AVX equivalence contract.

### Independent pairing-friendly-curves vectors

Retain canonical:

- G1 and G2 generator coordinates;
- compressed generator encodings;
- twelve Fp coefficients for the generator pairing result; and
- serialization mapping `[0,1,4,5,8,9,2,3,6,7,10,11]`.

Record the exact draft version and compare library final exponentiation
directly, avoiding the Python harness's generic exponentiation in routine tests.

## Planned-removal and CI triage

### Closed as literal removal-only

| Commits | Reason |
| --- | --- |
| `3aa4b71`, `374b7cd`, `097818a`, `ce49aad`, `7c17bf5`, `3ede151` | Go documentation, lint, WASM stub, or workflow only; removed in Phase 6B |
| `b79e2bc`, `2594ed1` | Python/SWIG-only compatibility; removed in Phase 6C |
| `d497208`, `0438d74` | Java Android change and its revert |
| `7f04257`, `68c1d90`, `8065152`, `4eb6250`, `f41f715`, `8aaaef0` | Zig binding/build/generator/CI/version only |
| `dafa98f` | Emscripten workflow and Zig log text only; local Emscripten JavaScript surface removed in Phase 6C |
| `c1bb460`, `6056ba3`, `9a5af53`, `d29fef4` | TypeScript/Node/SWIG Node-API only; removed in Phase 6C |

These rows are `not_applicable`, not deleted from the ledger. Their files
remain in the file manifest until reviewed removal is committed.

### `17226b2`: compressed infinity test intent

Decision: `translate`.

The upstream Go production code did not change; the test previously used the
wrong encoded length and therefore did not exercise the intended G2
uncompress path. Add direct retained tests for:

- G1 compressed infinity with 48 bytes;
- G2 compressed infinity with 96 bytes;
- scalar C/C++, Rust, and AVX public wrappers.

AVX wrappers are relevant because decoding delegates to scalar code.

### `85808fd` and `33d3f3a`: Rust WASI CI

Decision: `translate` as one retained workflow task.

If CatBLST advertises Rust WASI compatibility, test `wasm32-wasip1`. Cache a
pinned, checksum/signature-verified Wasmtime artifact. Do not preserve the
upstream unpinned installation pattern.

### `aee68f7`: empty salt

Decision: `not_applicable` for Go, closed with a retained link.

Rust empty slices already provide the non-null/zero-length case that triggers
`AUD-004`. Its C and Rust regression belongs to `a3c15d6`; Go code was removed
in Phase 6B.

### `06a5711`: retained uniqueness tree

Decision: `translate` as an audit/fix task.

`blst/common/rb_tree.cpp` uses fixed arrays sized to pointer bits for traversal
nodes and directions, but the loop has no explicit bound before indexing.
Do not blindly copy Go's return-false behavior because it conflates duplicate
input with resource exhaustion.

Audit:

- provable red-black maximum depth on 32- and 64-bit targets;
- `blst_uniq_sizeof` multiplication/addition overflow;
- allocation alignment and Rust byte-to-word sizing;
- large node counts and duplicate behavior; and
- a distinguishable failure contract if a limit is reachable.

### `f62244e`: retained workflow maintenance

Decision: `translate`.

Update checkout/cache/CodeQL actions during CI hardening using deliberately
selected current versions. Drop Go/Node portions and verify downloaded
toolchains. Do not treat action-version refreshes as security coverage by
themselves.

### `db3defd`: corroborating Pippenger evidence

Decision: `not_applicable` for Go, closed with a retained link.

The Go fix independently confirms the small-input/high-worker divide-by-zero
class already recorded under Rust `6271918` and `AUD-003`. Required retained
tests include windows 2/3, 0/1 points, points fewer than workers, large mocked
worker counts, and both G1/G2 paths.

## Wave 02 result

Every upstream-only commit now has a non-pending decision. Rows remain
`reviewing` when implementation, policy, or validation is outstanding. Rows are
`closed` only when the literal change is proven removal-only, deliberately
rejected with rationale, or already equivalent with validation.
