# File Review: SHA-256 streaming header

## Record

| Field | Value |
| --- | --- |
| Path | `blst/sign/sha256.h` |
| Roadmap phase / queue | Phase 15 scalar curves and protocols, eighth packet |
| Reviewed commit | `5c69ad63c76d9e63e9cf75a708ac84ede45b5eb6` |
| Reviewer | Codex audit team |
| Date | 2026-08-05 |
| Disposition | in-scope production |
| Subsystem | SHA-256 one-shot and streaming state used by keygen and hash-to-field |
| Profiles | scalar selected SHA backend; C/C++/generated Rust one-shot boundary |
| Generation/provenance | CatBLST C++ adaptation of pinned BLST `src/sha256.h` |
| Related files | `blst/exports.cpp`, `blst/sign/keygen.cpp`, `blst/scalar/fields/hash_to_field.h`, `bindings/blst_scalar.h`, `test/Sha256_audit_tests.cpp` |

## Purpose and trust boundary

The header owns SHA-256 state initialization, streaming updates, padding, and
digest publication. Messages can contain public protocol data or secret keygen
material. The fixed compression backend is selected at compile time. The native
streaming API is internal C++; `blst_sha256` is the retained public allocation-
free C one-shot boundary, and generated raw Rust inherits its C contract.

## Contracts and invariants

- `Sha256Context` is initialized before update/final and reinitialized before
  reuse after final. Direct field modification is unsupported.
- Input pointers cover their stated byte ranges and may be null only for zero
  length. Digest pointers cover 32 writable bytes.
- The supported message domain is fewer than `2^64` bits. Debug verifies the
  accumulated count, buffered count, and their modulo-64 relation.
- Streaming operations are allocation-free and `noexcept`. Their control flow
  depends on public message length and buffered count, not message contents.
- Internal finalization consumes but does not erase caller-owned state. Owners
  explicitly call `sha256_clear`; the public one-shot wrapper always does so
  after digest publication.
- Public one-shot output may overlap input because all input is consumed before
  digest emission.
- Functional zero filling is ordinary optimization-friendly storage. Secure
  lifecycle erasure remains the volatile `vec_zero` primitive.
- `SHA256_CTX` remains only as an internal compatibility alias for frozen
  BTx/BTE; new supported code uses `Sha256Context`.

## Dependency and call-path map

`blst_sha256` creates, updates, finalizes, and erases one local context. Keygen
uses the state for HMAC, HKDF, and EIP-2333, and erases its complete owning
scratch. Hash-to-field uses the same streaming state plus the standard XMD
precomputed Z-pad state. BTx/BTE includes the alias but remains frozen and was
neither edited nor audited. Backend compression, emit, and copy primitives were
reviewed previously under `AUD-055`.

## Upstream/specification comparison

The initial hash constants, 64-byte block schedule, padding cutoff, big-endian
bit-length encoding, and compression sequence match pinned BLST. CatBLST adds
fixed-width/state assertions, explicit lifecycle documentation, a supported
public declaration, independent end-to-end tests, fast functional zeroing, and
one-shot context erasure. No SHA arithmetic or backend dispatch was changed.

## Checklist evidence

1. Every declaration, state field, branch, direct caller, public export, and
   generated raw Rust declaration was traced.
2. Pointer, length, overlap, lifecycle, allocation, exception, timing, and
   erasure contracts are explicit.
3. NIST empty, `abc`, multi-block, and million-`a` vectors pass.
4. Independent 55/56/63/64/65/119/120/127/128-byte vectors pass one-shot,
   byte-at-a-time, and mixed streaming schedules.
5. Null-empty input, exact/partial output overlap, explicit full-context clear,
   and Debug invariant failures are retained regressions.
6. Existing keygen/EIP-2333 and XMD/map regressions pass 13/13.
7. MSVC C11 `/W4 /WX` compiles the public declaration; the C++ target links it.
8. MSVC 19.50 optimized output uses eight 64-bit functional-zero stores and a
   fixed 112-byte secure one-shot erase after digest emission.
9. The optimized object contains one bare public SHA symbol and no bare internal
   SHA helper symbols.
10. Audit manifest and public-documentation checkers pass.

## Findings and disposition

`AUD-141` through `AUD-145` are validated. No critical/high residual defect is
open in this file.

## Residual risk and deferred evidence

Sanitizer/fuzz coverage, other host compiler code generation, Apple execution,
and the full supported-language milestone matrix remain release evidence. The
experimental CUDA path is outside this gate. Compiler-created copies outside
the explicit context object remain beyond the C++ abstract-machine erasure
guarantee. HMAC/HKDF/keygen semantics and complete owning-scratch lifecycle are
reviewed with `blst/sign/keygen.cpp` next.

## Completion

Complete Phase 15 file review at implementation commit
`5c69ad63c76d9e63e9cf75a708ac84ede45b5eb6`.

## Batch E backend follow-up

The streaming interface and padding logic are unchanged at `2f59883`. Its
selected compression boundary now has an explicit provider chain: SHA-NI,
SSSE3, then portable C on x86; hardware SHA then NEON on ARM64. Every provider
retains the same state and block contract. Focused known-answer tests exercise
each executable x86 provider. Supported Apple instructions are unchanged from
the prior hosted run.
