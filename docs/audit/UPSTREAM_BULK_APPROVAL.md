# Upstream Bulk Approval

## Record

- Date: 2026-07-31
- Scope: the final 10 upstream rows originally marked `reviewing`
- Decision state: approved for implementation under the translations below
- Validation state: `07b109d` is superseded and closed; the three CHERI rows
  are validated; six rows remain
- Benchmarks: excluded and not run
- BTx/BTE: frozen and untouched
- CUDA: retained experimental and outside these packets

This record completes the decision gate for the upstream manifest. Approval is
not validation: outstanding rows remain `reviewing` until their implementation
packets have focused evidence and the manifest can truthfully move to
`validated` or `closed`.

## Approved implementation packets

| Packet | Rows | Approved translation | Required focused evidence |
| --- | --- | --- | --- |
| Caller-threaded Rust runtime | `07b109d` | Superseded after architectural review: reject the bounded-channel patch and remove all internal Rust scheduling. Preserve the optimized native arithmetic kernels while every public call executes entirely on its caller thread. | Default/full/portable/no-std Rust evidence; deterministic errors; dependency and source-policy proof that no internal scheduler remains; see `PHASE_08.md`. |
| CHERI inverse scratch bounds | `2e1948b`, `ed865dc`, `fe85dc3` | Port `scbnds` support and both inverse-generator changes, then regenerate the exact 10 affected ARM artifacts through the verified assembly pipeline. Do not import the incidental `blst.def` delta. | Deterministic regeneration and exact allowlist; ordinary ARM64 assembly/compile and inverse tests; CHERI assembly/execution when available, otherwise a recorded hardware/toolchain block. |
| Retained CI security and WASI | `85808fd`, `33d3f3a`, `68b58bd`, `f62244e` | Implement one current-policy workflow update: retained Rust `wasm32-wasip1` coverage using a versioned, checksum-verified Wasmtime artifact; CodeQL v4 for C/C++ and Rust; current retained actions pinned to full commit SHAs. Drop all Go, Node, SWIG, C#, and Zig portions. | Workflow syntax; runner compatibility; C/C++ manual CodeQL build; Rust CodeQL coverage without unsupported bindings; WASI build and execution; action/download provenance. |
| Conformance vectors and README language | `ff66582`, `67fc54b` | Move the draft-13 BLS12-381 generator coordinates, compressed encodings, and pairing coefficients into retained C++ and Rust tests; exercise scalar and AVX-512 production paths where those paths exist; do not restore Python. Translate the README intent, not its literal patch: do not add Zig, and replace the floating V5 statement only after retained-language signature vectors establish the exact supported revision. | Exact pairing-draft provenance and serialization mapping; explicit treatment of BLST's production pairing cofactor convention; scalar/AVX differential coverage; C++ and Rust signature vectors for both group orientations; an exact version-pinned Internet-Draft claim matching the evidence. |

## Rationale and exclusions

### Rust scheduling amendment

The original bulk review approved a partial bounded-channel translation.
Architectural review then established that the process-global CPU-sized pool,
unsafe scoped-job lifetime erasure, and hidden scheduling contract were not
acceptable for CatBLST. Phase 8 removes every internal queue and scheduler, so
`07b109d` is now deliberately rejected as obsolete rather than implemented.

### CHERI generated outputs

The upstream generated commit also contains two Win64 exports associated with
an API change already translated in Phase 7A. `AUD-016` reserves the complete
Windows export surface for ABI/prefix hardening, so this packet must not
partially patch `build/win64/blst.def`.

### CI supply chain

The upstream Wasmtime installation used an unversioned network installer. The
retained job must instead consume a specific release artifact and verify its
digest. GitHub Actions references must be immutable full SHAs, with the human
readable release version recorded alongside them. CodeQL is supplementary
static analysis; it does not cover handwritten assembly, constant-time
properties, or scalar/AVX equivalence.

This approval explicitly retains `wasm32-wasip1` as a tested target of the
supported Rust crate. It does not restore the removed Emscripten/JavaScript
binding. The crate already advertises WebAssembly compilation and is now
caller-threaded on every target, so an executable CI target is the appropriate
way to keep that claim honest.

### Draft scope

The `ff66582` vectors come from
[`draft-irtf-cfrg-pairing-friendly-curves-13`](https://datatracker.ietf.org/doc/draft-irtf-cfrg-pairing-friendly-curves/13/),
an active Internet-Draft dated 2026-07-06. They are valuable independent curve,
serialization, and pairing tests, but they do not establish BLS signature
draft conformance by themselves. The draft also documents that production
pairing libraries may use a cofactor convention whose raw output differs from
the literal pseudocode, so the test must state which value it compares.

The literal `67fc54b` patch changed V5 to V6 and added Zig. Zig is absent. As of
this review, [`draft-irtf-cfrg-bls-signature-06`](https://datatracker.ietf.org/doc/draft-irtf-cfrg-bls-signature/06/)
is the latest published revision but expired on 2026-05-06. A versioned
historical compatibility claim may still be useful, but only after its own
signature vectors pass; the README must not imply that an expired
Internet-Draft is an IETF standard.

## Execution order

1. Caller-threaded Rust runtime - completed in Phase 8.
2. CHERI/ARM generator and exact regenerated artifacts - completed in Phase 9.
3. Pairing and signature conformance evidence, followed by accurate README
   language.
4. CI/WASI/CodeQL policy update, after the test targets it should invoke are
   stable.

The packets are bulk-approved in this order and do not require another design
approval unless implementation uncovers a material scope or API change.
