# File Review: `blst/asm/ctx_inverse_mod_384-x86_64.pl`

## Record

| Field | Value |
| --- | --- |
| Reviewed commit | `83c11ba829a018e1ce38297601ac40286f2d6791` |
| Queue | Phase 13 assembly Batch D |
| Status | `complete` |
| Upstream | `54e6e55674722fc2797ebb4bbb71b26d881eb4b8` |
| Findings | `AUD-372`, `AUD-373` |

## Review

This owner emits the ADX 384-bit almost-inverse. Its arithmetic matches pinned
upstream. Twenty-three 31-step groups and a 55-step tail perform 768 divsteps.
All input choices use masks and conditional moves. Branches count fixed loops.
The 1,112-byte frame contains one aligned 512-byte scratch region.

The direct ADX independent byte oracle passes the exact double-width bound and
congruence for zero, boundary, and dense inputs. Generated COFF, ELF, Mach-O,
and Win64 outputs reproduce deterministically. No benchmark was required.
