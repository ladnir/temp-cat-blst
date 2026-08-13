# File Review: `blst/asm/ctq_inverse_mod_384-x86_64.pl`

## Record

| Field | Value |
| --- | --- |
| Reviewed commit | `83c11ba829a018e1ce38297601ac40286f2d6791` |
| Queue | Phase 13 assembly Batch D |
| Status | `complete` |
| Upstream | `54e6e55674722fc2797ebb4bbb71b26d881eb4b8` |
| Findings | `AUD-372`, `AUD-373` |

## Review

This owner emits the generic 384-bit x86-64 provider and publicly dispatches
to the ADX provider when the public capability bit is set. The only local
semantic difference from pinned upstream is the approved capability-symbol
spelling. The generic body performs twelve 62-step groups and a 24-step tail.
Secret choices are branchless; the provider branch and loop counters are
public.

The forced-generic independent byte oracle validates the result bound and
`T * inp = 2^768 (mod p)` relation. Deterministic generation reproduces every
x86 dialect. No instruction or benchmark changed.
