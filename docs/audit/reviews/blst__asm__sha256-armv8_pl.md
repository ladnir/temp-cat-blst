# File Review: `blst/asm/sha256-armv8.pl`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/asm/sha256-armv8.pl` |
| Roadmap phase / queue | Phase 13, scalar assembly trust chain, Batch E |
| Reviewed commit | `2f59883caf67191aae95abb905598f405f8139e8` |
| Reviewer | Codex audit team |
| Date | 2026-08-13 |
| Disposition | `in_scope_production` |
| Subsystem | AArch64 SHA-256 compression and provider dispatch |
| Profiles | supported Apple ARM64; retained experimental ARM64 and CHERI |
| Upstream reference | `54e6e55674722fc2797ebb4bbb71b26d881eb4b8:src/asm/sha256-armv8.pl` |
| Related files | `common/cpuid.cpp`, `backend_asm.h`, five generated ARM outputs, `test/Sha256_audit_tests.cpp` |

## Purpose and contracts

The owner generates one public compression entry with two implementations.
The public capability word selects ARMv8 SHA instructions when available and
the generic NEON schedule otherwise. Both implementations accept a positive
public block count, an eight-word state, and consecutive 64-byte blocks. State
and message storage are disjoint. Message contents do not select branches or
addresses.

Apple ARM64 is supported. Other AArch64 dialects and CHERI remain experimental
and do not contribute release evidence.

## Review evidence

The hardware rounds, generic NEON schedule, endian conversion, state
accumulation, block traversal, stack frame, ABI translations, and helper
functions were read against the upstream owner. No SHA arithmetic changed.
The supported Apple initializer selects the hardware provider. Retained hosts
without the hardware capability use the NEON fallback.

`AUD-376` found that the CHERI generic schedule recovered pointers from the
stack capability without restricting their authority to the 64-byte schedule
region. The owner now bounds each recovered schedule pointer. All current
accesses were already in range, so this is confinement hardening rather than a
confirmed wrong-result repair.

The CHERI translator fixture passes. Deterministic generation reproduces all
97 outputs. Aggregate ARM64 ELF, COFF, and Mach-O object translation passes.
Preprocessing removes every CHERI-only bound from supported Apple code. The
nonblank Apple instruction projection is identical to the reviewed baseline.
Prior supported Apple Silicon SHA known-answer execution therefore remains
applicable.

## Findings and decision

- `AUD-376` — accepted risk: CHERI schedule capabilities are restricted to 64
  bytes, but native execution remains mandatory before experimental-target
  promotion.
- `AUD-377` — validated: block, overlap, and timing contracts are explicit.
- `AUD-378` — validated for supported Apple ARM64 through unchanged projection
  and prior hosted known-answer execution.

The supported Apple review is complete at the recorded commit. CHERI remains
experimental until native assembly and execution validate the bounded
capabilities. No performance change reaches Apple ARM64.
