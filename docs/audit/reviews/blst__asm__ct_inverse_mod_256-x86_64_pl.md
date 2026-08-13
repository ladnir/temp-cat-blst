# File Review: `blst/asm/ct_inverse_mod_256-x86_64.pl`

## Record

| Field | Value |
| --- | --- |
| Reviewed commit | `83c11ba829a018e1ce38297601ac40286f2d6791` |
| Queue | Phase 13 assembly Batch D |
| Status | `complete` |
| Upstream | `54e6e55674722fc2797ebb4bbb71b26d881eb4b8` |
| Findings | `AUD-372`, `AUD-373` |

## Review

This owner emits the 256-bit x86-64 fixed-schedule almost-inverse. Its
arithmetic is unchanged from pinned upstream. Local differences are the
audited argument-list translator invocation and failure propagation.

The routine performs fifteen 31-step groups and one 47-step tail, for 512
divsteps. Input-dependent choices use arithmetic masks. Branches count fixed
loops. The 1,072-byte frame contains one aligned 512-byte scratch region and
fixed spill storage. The result is a 512-bit value `T` with
`T * inp = 2^512 (mod mod)` for invertible nonzero input. Zero returns zero.

Selected x86 execution and the independent byte oracle pass zero, one,
`r-1`, and a dense reduced value. Deterministic generation reproduces all four
x86 dialects. No instruction or benchmark changed.
