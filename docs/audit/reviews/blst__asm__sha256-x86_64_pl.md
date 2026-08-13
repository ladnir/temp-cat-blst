# File Review: `blst/asm/sha256-x86_64.pl`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/asm/sha256-x86_64.pl` |
| Roadmap phase / queue | Phase 13, scalar assembly trust chain, Batch E |
| Reviewed commit | `2f59883caf67191aae95abb905598f405f8139e8` |
| Reviewer | Codex audit team |
| Date | 2026-08-13 |
| Disposition | `in_scope_production` |
| Subsystem | x86-64 SHA-256 compression and provider dispatch |
| Profiles | x86-64 ELF, COFF, Mach-O, and Win64 |
| Upstream reference | `54e6e55674722fc2797ebb4bbb71b26d881eb4b8:src/asm/sha256-x86_64.pl` |
| Related files | `common/cpuid.cpp`, `backend_portable.cpp`, `backend_asm.h`, generated x86 outputs, `test/Sha256_audit_tests.cpp` |

## Purpose and contracts

The owner generates SHA-NI and SSSE3 compression providers plus one public
dispatcher. Compression accepts a positive public block count, an eight-word
state, and consecutive 64-byte blocks. State and message storage are disjoint.
The routine allocates nothing and does not throw. Control flow and addresses
depend on the block count and public CPU capabilities, not message contents.

The selected provider order is SHA-NI, SSSE3, then portable C. The dispatcher
uses direct tail jumps and creates no stack frame. SGX-LVI builds retain the
inherited restriction that disables SHA-NI selection.

## Review evidence

The SHA-NI and SSSE3 round schedules, endian conversion, state accumulation,
block traversal, ABI prologues, and returns were read against the upstream
owner. Their arithmetic is unchanged. `AUD-375` records the prior unsupported
assumption that every non-SHA-NI x86 CPU had SSSE3. The unreachable scalar
x86 generator did not repair that assumption.

CPUID now records SSSE3 separately. The public entry tests SHA-NI once, then
SSSE3, and otherwise tail-jumps to the portable C compression routine. The
optimized Win64 object contains the expected stackless selector and no
conversion or out-of-line wrapper before the provider entry.

Forced provider tests use the same `abc` and deterministic two-block vectors
for C, SSSE3, and SHA-NI. SHA-NI is forced only when the executing host reports
the instruction set. Six focused SHA tests pass under MSVC 19.50 Release. Five
pass in the nonportable ADX configuration. Perl syntax, deterministic
generation of 97 outputs, and aggregate translator-object checks pass.

## Findings and decision

- `AUD-375` — validated: every x86 capability state has a defined provider.
- `AUD-377` — validated: block, overlap, and timing contracts are explicit.
- `AUD-378` — validated: each executable x86 provider has direct KAT evidence.

The review is complete at the recorded commit. Future x86 providers must add a
distinct capability bit and direct known-answer coverage. Benchmarking was not
required because the preferred SHA-NI body is unchanged and gains only one
public capability test plus a direct tail jump.
