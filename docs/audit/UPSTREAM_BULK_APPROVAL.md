# Upstream Bulk Approval

## Record

- Date: 2026-07-31
- Scope: the final 10 upstream rows originally marked `reviewing`
- Decision state: approved for implementation under the translations below
- Validation state: the runtime, CHERI, and conformance packets are complete;
  WASI was removed from supported scope; two retained CI/security rows remain
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
| Retained CI security and provenance | `68b58bd`, `f62244e` | Implement one current-policy workflow update: CodeQL v4 for C/C++ and Rust, plus current retained actions pinned to full commit SHAs. Drop all Go, Node, SWIG, C#, Zig, and WASI portions. | Workflow syntax; runner compatibility; C/C++ manual CodeQL build; Rust CodeQL coverage without unsupported bindings; action provenance; integration with the standing automated-analysis workstream in `AUTOMATED_ANALYSIS_PLAN.md` and `PRODUCTION_AUDIT_ROADMAP.md`. |
| Conformance vectors and README language | `ff66582`, `67fc54b` | Move the draft-13 BLS12-381 generator coordinates, compressed encodings, and pairing coefficients into retained C++ and Rust tests; exercise scalar and AVX-512 production paths where those paths exist; do not restore Python. Translate the README intent, not its literal patch: do not add Zig, and replace the floating V5 statement only after retained-language signature vectors establish the exact supported revision. | Completed in Phase 11: exact pairing-draft provenance and serialization mapping; explicit cubed production-pairing convention; scalar/AVX coverage; C++ and Rust signature vectors for both group orientations; narrow version-pinned draft-07 language. |

## Rationale and exclusions

### Rust scheduling amendment

The original bulk review approved a partial bounded-channel translation.
Architectural review then established that the process-global CPU-sized pool,
unsafe scoped-job lifetime erasure, and hidden scheduling contract were not
acceptable for CatBLST. Phase 8 removes every internal queue and scheduler, so
`07b109d` is now deliberately rejected as obsolete rather than implemented.

### CHERI generated outputs

The upstream generated commit also contains two Win64 exports associated with
an API change already translated in Phase 7A. At this packet, `AUD-016`
reserved the complete Windows export surface, so the packet did not partially
patch `build/win64/blst.def`. The later closure regenerated the full surface.

### CI supply chain and WASI amendment

GitHub Actions references must be immutable full SHAs, with the human-readable
release version recorded alongside them. CodeQL is supplementary static
analysis; it does not cover handwritten assembly, constant-time properties, or
scalar/AVX equivalence.

Architectural review removed `wasm32-wasip1` from the supported production
matrix. CatBLST is focused on native C, C++, and Rust over scalar and AVX-512,
with Apple Silicon in supported scalar scope. The Wasmtime cache and its YAML
repair therefore add an unsupported platform and supply-chain surface without
advancing the primary release gate. Rows `85808fd` and `33d3f3a` are closed as
rejected/not applicable; any remaining WebAssembly compatibility language must
be removed or narrowed during the relevant documentation review.

### Draft scope

The `ff66582` vectors come from
[`draft-irtf-cfrg-pairing-friendly-curves-13`](https://datatracker.ietf.org/doc/draft-irtf-cfrg-pairing-friendly-curves/13/),
an active Internet-Draft dated 2026-07-06. They are valuable independent curve,
serialization, and pairing tests, but they do not establish BLS signature
draft conformance by themselves. The draft also documents that production
pairing libraries may use a cofactor convention whose raw output differs from
the literal pseudocode, so the test must state which value it compares.

The literal `67fc54b` patch changed V5 to V6 and added Zig. Zig is absent.
[`draft-irtf-cfrg-bls-signature-07`](https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-bls-signature-07)
superseded draft-06 on 2026-07-06. Because its test-vector appendix remains
unfinished, Phase 11 cross-checked the two basic-scheme group orientations
against pinned `@noble/curves` 2.2.0 output. The README now makes only that
narrow, versioned ciphersuite statement and explicitly identifies the document
as a work-in-progress Internet-Draft rather than an IETF standard.

## Execution order

1. Caller-threaded Rust runtime - completed in Phase 8.
2. CHERI/ARM generator and exact regenerated artifacts - completed in Phase 9.
3. Pairing and signature conformance evidence, followed by accurate README
   language - completed in Phase 11.
4. Initial audit-platform milestone followed by the standing automated-analysis
   workstream and retained CodeQL/action provenance; WASI is excluded.

The packets are bulk-approved in this order and do not require another design
approval unless implementation uncovers a material scope or API change.
