# File Review: `blst/asm/ct_inverse_mod_384-armv8.pl`

## Record

| Field | Value |
| --- | --- |
| Reviewed commit | `83c11ba829a018e1ce38297601ac40286f2d6791` |
| Queue | Phase 13 assembly Batch D |
| Status | `complete` |
| Upstream | `54e6e55674722fc2797ebb4bbb71b26d881eb4b8` |
| Findings | `AUD-372`, `AUD-373` |

## Review

This owner emits the generic AArch64 384-bit almost-inverse. Its arithmetic
matches pinned upstream. Twelve 62-step groups and a 24-step tail perform 768
divsteps. Conditional selects implement data choices; loop branches have fixed
counts. The fixed 1,056-byte frame contains one aligned, CHERI-bounded
512-byte scratch region.

Supported Apple Silicon scalar-assembly job `94666845914` in run `31767696988`
compiles this owner and passes the complete CTest graph, including the
independent double-width byte oracle. The supported Apple review is complete.
Non-Apple ARM64 and CHERI remain experimental.
