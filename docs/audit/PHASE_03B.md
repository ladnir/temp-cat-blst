# Phase 3B Win64 Assembly Symbol Evidence

## Record

- Date: 2026-07-30
- Baseline commit: `a204b4c`
- Implementation commit: `0bceafe`
- State: validated and committed
- Scope: `AUD-010`, upstream `41b511a`, and generated companion `14cd09a`
- Benchmarks: disabled and not run
- CUDA: disabled

## Root cause

The x86 generators emitted Win64 `EXTERN` declarations for q-to-x fallback
symbols and matching `PUBLIC` declarations for the `$1` provider entry points
unconditionally. The runtime branches were guarded, but the object symbol
tables were not. Non-portable MSVC builds therefore exposed 20 unnecessary
fallback providers and retained 20 matching undefined consumer references.

The repository also used two spellings for the same assembly policy:
`__BLST_PORTABLE__` in the generators and `BLST_PORTABLE` in two manually
divergent generated MASM files. CMake did not pass either spelling to MASM.

## Policy and implementation

`BLST_PORTABLE` remains the public CMake option. CMake now passes
`__BLST_PORTABLE__` explicitly to every MASM source when that option is
enabled, matching the existing translation in `build/assembly.S`.

The six affected x86 generators now guard both sides of the Win64 fallback
linkage:

- q consumers publish `$1` references only in portable mode;
- x providers publish the corresponding `$1` definitions only in portable
  mode.

The temporary Phase 3A MASM-only spelling substitutions were removed.
Regeneration changed only the six expected Win64 outputs. This is the local
translation of upstream `41b511a` plus its generated-output companion
`14cd09a`; it does not yet import the other accepted assembly schedule or
translator changes assigned to Phase 3C.

## Symbol allowlist

`scripts/audit/check_win64_asm_symbols.ps1` inspects every MASM object with
`dumpbin /symbols` and requires an exact mode-specific surface.

Portable mode contains exactly 40 external fallback entries:

- 20 undefined entries in the three q consumer objects;
- 20 defined entries in the three x provider objects.

Non-portable mode contains no external symbol ending in `$1` in any MASM
object. Unexpected symbols, missing symbols, incorrect provider/consumer
placement, or incorrect definition state fail the check.

## Validation

Toolchain:

- CMake with Ninja;
- MSVC C/C++ 19.50.35728.0;
- MASM 14.50.35728.0;
- `dumpbin` from MSVC 14.50.35717.

| Profile | Configuration | Symbol result | Test result |
| --- | --- | --- | --- |
| Portable dispatch | Debug, scalar default, `BLST_PORTABLE=ON`, `BLST_ADX=OFF`, assembly on | Exact 40-entry allowlist across 13 MASM objects | 78/78 focused interface, field, Montgomery, and shift tests passed |
| Non-portable ADX | Debug, scalar default, `BLST_PORTABLE=OFF`, `BLST_ADX=ON`, assembly on | Exact zero-entry allowlist across 13 MASM objects | 38/38 focused field, Montgomery, shift, group, and pairing tests passed |

`build/verify-assembly.sh` generated two independent trees and verified all
100 outputs against the checked-in allowlist. `git diff --check` passed.

The portable profile deliberately disabled compile-time ADX so calls entered
the q functions and exercised their CPU-dispatch path. The non-portable
profile enabled ADX for direct optimized entry points.

## Residual validation

- Validation covered Win64 x86-64, the affected platform and architecture.
- Release-mode code generation was not separately tested; the guarded MASM
  symbol surface is configuration-independent.
- The complete test suite was not rerun because Phase 3A already ran the full
  scalar and AVX-portable suites. Phase 3B ran the assembly-sensitive focused
  sets in both modes.
- The remaining accepted upstream assembly generator, schedule, comment, and
  translator changes remain Phase 3C work.
