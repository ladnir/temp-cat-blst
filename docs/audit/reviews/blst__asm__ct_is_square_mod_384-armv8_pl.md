# File Review: `blst/asm/ct_is_square_mod_384-armv8.pl`

## Record

| Field | Value |
| --- | --- |
| Reviewed commit | `83c11ba829a018e1ce38297601ac40286f2d6791` |
| Queue | Phase 13 assembly Batch D |
| Status | `complete` |
| Upstream | `54e6e55674722fc2797ebb4bbb71b26d881eb4b8` |
| Findings | `AUD-372`, `AUD-373`, `AUD-374` |

## Review

This owner is the AArch64 residue test. Its arithmetic matches the x86 owner
and pinned upstream: 24 fixed 30-step groups plus a 48-step tail. Conditional
selects implement data choices, and loop branches have fixed counts.

`AUD-374` bounds the aligned CHERI scratch capability to exactly 256 bytes.
Each alternating pointer derives from a bounded peer, so XOR bit 7 cannot
restore authority over the full stack. All five outputs regenerate exactly.
The ordinary Apple preprocessed instruction stream is unchanged. Supported
Apple Silicon scalar-assembly job `94666845914` in run `31767696988` compiles
this owner and passes the independent Euler-oracle regression. The supported
Apple review is complete. CHERI stays experimental pending Morello assembly
and execution.
