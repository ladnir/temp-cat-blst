# Phase 7A Variable-Length Fp Decoder Evidence

## Record

- Date: 2026-07-31
- Implementation commit: `f962830`
- State: validated and committed
- Scope: upstream `cf739f0` and generated companion `9c35e86`
- Benchmarks: excluded and not run
- CUDA: experimental and not exercised
- BTx/BTE: frozen and untouched

## Decision

Translate the upstream arbitrary-length Fp decoding API into the fork's split
C header and native C++ export architecture. The implementation remains
allocation-free and processes 48-byte radix chunks with one Montgomery
multiplication per chunk.

The local contract is deliberately stronger for empty input: both functions
set the output to field zero and accept `in == NULL` when `len == 0`. This
avoids the null-pointer arithmetic present in the upstream little-endian empty
case. A nonzero length still requires a readable input buffer.

The Rust FFI declarations were updated with the exact generated companion
surface. The bindgen CLI was not available, so this six-line companion was
reviewed and applied directly rather than claiming a full regeneration.

## Validation

The focused validation required by `VALIDATION_POLICY.md` passed:

| Surface | Focused evidence | Result |
| --- | --- | --- |
| C/native C++ direct scalar | New bytewise-oracle regression covering null-empty input, chunk boundaries, long inputs, reversal, and all-ones inputs | Passed |
| Portable scalar dispatch | Same focused decoder regression in the portable-dispatch build | Passed |
| Rust FFI | Focused empty-input and 49-byte endian-equivalence/link test | Passed |
| C and C++ headers | Standalone C11 and C++23 syntax compilation of both declarations | Passed |
| Exported ABI | MSVC archive symbol inspection | Both `blst_fp_from_le_bytes` and `blst_fp_from_be_bytes` are external |
| AVX-512 | Existing AVX-512 target rebuilt | Passed; no AVX implementation changed |
| Formatting/integrity | `cargo fmt --check` and `git diff --check` | Passed |

Before the proportional policy was adopted, the direct scalar,
portable-dispatch, and Rust full suites also passed with zero failures. A
duplicate full AVX-512 sweep was stopped because it added no useful evidence
for this scalar export change. No sanitizer sweep was run; focused boundary
coverage and the explicit zero-length guard exercise the changed pointer
contract, while sanitizers remain part of the next relevant subsystem or
milestone gate.

## Residual risk

- At this phase, `build/win64/blst.def` remained globally stale under
  `AUD-016`; the later export-hardening packet regenerated the complete
  compatibility surface and added a deterministic CI comparison.
- The generated Rust declaration update should be reproduced when the retained
  bindgen workflow is next exercised.
- No ARM execution was performed; the implementation uses the existing
  portable/assembly backend selection already exercised by scalar builds.
