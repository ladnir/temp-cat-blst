# File Review: `blst/asm/ct_is_square_mod_384-x86_64.pl`

## Record

| Field | Value |
| --- | --- |
| Reviewed commit | `83c11ba829a018e1ce38297601ac40286f2d6791` |
| Queue | Phase 13 assembly Batch D |
| Status | `complete` |
| Upstream | `54e6e55674722fc2797ebb4bbb71b26d881eb4b8` |
| Findings | `AUD-372`, `AUD-373` |

## Review

This owner computes quadratic residuosity modulo an odd prime. Twenty-four
30-step groups and one 48-step tail perform 768 divsteps. Input choices use
masks and conditional moves. The only branches count fixed loops. The fixed
536-byte frame contains one aligned 256-byte scratch region.

An independent Euler-criterion oracle confirms that zero and residues return
true and nonresidues return false. The portable result is checked separately.
All four x86 output dialects reproduce deterministically, and no instruction
or benchmark changed.
