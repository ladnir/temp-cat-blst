# File Review: `blst/asm/ct_inverse_mod_256-armv8.pl`

## Record

| Field | Value |
| --- | --- |
| Reviewed commit | `83c11ba829a018e1ce38297601ac40286f2d6791` |
| Queue | Phase 13 assembly Batch D |
| Status | `complete` |
| Upstream | `54e6e55674722fc2797ebb4bbb71b26d881eb4b8` |
| Findings | `AUD-372`, `AUD-373` |

## Review

This owner is the AArch64 counterpart of the 256-bit almost-inverse. The
arithmetic matches pinned upstream. It performs the same fixed 512 divsteps,
uses conditional selects instead of secret branches, and reserves a fixed
1,040-byte frame.

The aligned 512-byte scratch capability is bounded on CHERI. XOR bit 8 selects
its two halves, and derived capabilities preserve the bound. Supported Apple
Silicon scalar-assembly job `94666845914` in run `31767696988` compiles this
owner and passes the complete CTest graph, including the independent byte
oracle. The supported Apple review is complete. Non-Apple ARM64 and CHERI stay
experimental.
