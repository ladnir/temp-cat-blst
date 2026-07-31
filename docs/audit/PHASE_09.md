# Phase 9 CHERI Inverse Scratch Evidence

## Record

- Date: 2026-07-31
- Implementation commit: `ce8e666`
- State: validated with explicit target-toolchain residuals
- Scope: upstream `2e1948b`, `ed865dc`, `fe85dc3`; `AUD-022`; and
  `AUD-023`
- Benchmarks: not run or audited
- CUDA: experimental and not exercised
- BTx/BTE: frozen and untouched

## Decision and implementation

The three upstream changes were translated through the verified local
generator path rather than cherry-picked across the fork's relocated source
layout.

- Added the CHERI `scbnds` synthetic instruction to `arm-xlate.pl`.
- Bounded the aligned inverse scratch capability to exactly 512 bytes in the
  256- and 384-bit ARM generators.
- Derived each flip-flop pointer from the already bounded peer capability
  instead of from `csp`, so changing its address cannot restore authority over
  the wider stack frame.
- Regenerated exactly the two inverse routines in the five checked-in ARM
  formats: CHERI, COFF, ELF, Mach-O, and Win64.
- Deliberately left `build/win64/blst.def` unchanged for `AUD-016`.

The 512-byte bound is sufficient and tight for the intended scratch region.
The 256-bit routine reserves 1,040 bytes and can place the aligned region at
frame offsets 16 through 512. The 384-bit routine reserves 1,056 bytes and can
place it at offsets 32 through 528. In both cases the complete 512-byte region
remains inside the frame. Its 512-byte alignment makes XOR bit 8 alternate
between its two 256-byte halves without leaving the bounded capability.

## Review findings fixed in the packet

`AUD-022` was found while reviewing the new translator handler. The third
operand conversion in `scvalue` and the new `scbnds` recognized only `c0`
through `c9`; `cadd` and `csub` handled larger registers only accidentally by
rewriting a single-digit prefix. All four now translate a complete register
number, and a permanent fixture covers x10 through x24, an immediate bound,
`csp`, CHERI output, and ordinary Linux output.

`AUD-023` was found while selecting the future ARM execution gate. The
256-bit inverse test passed the scalar modulus with a left-aligned base-field
modulus and compared only assembly with portable code. It now constructs the
correct scalar `modx`, checks zero, one, `r - 1`, and deterministic random
inputs, and verifies the multiplicative identity for both implementations.

## Focused validation

Following `VALIDATION_POLICY.md`, validation was limited to the affected
generator and inverse-test surfaces.

| Surface | Evidence | Result |
| --- | --- | --- |
| Perl sources | `perl -c` on the translator, both inverse generators, and the new fixture | Passed |
| Translator operands | `perl scripts/audit/check_arm_xlate_cheri.pl` | Passed for multi-digit CHERI operands, immediate bounds, `csp`, and non-CHERI `cadd`/`csub` |
| Assembly generation | `./build/refresh.sh verify-assembly` | All 100 outputs reproduced twice and matched the checked-in tree |
| Generated scope | Changed-path and upstream-blob comparison | Exactly 10 inverse artifacts changed; all 10 match upstream `fe85dc3` byte-for-byte |
| Ordinary ARM text | GNU preprocessing and Win64 false-branch projection against pre-change files | All 10 non-CHERI projections are byte-identical |
| Win64 exports | Direct Git comparison | `build/win64/blst.def` is unchanged |
| Inverse oracle | `InterfaceTests.CtInverseMod256` and `InterfaceTests.CtxInverseMod384` in the existing portable-plus-assembly scalar build | 2/2 passed |
| Formatting | `git diff --check` | Passed |

No full C/C++/Rust suite was warranted: the production delta is confined to
CHERI-guarded ARM assembly text, while the ordinary ARM projection is
unchanged. No benchmark was run.

## Residual target validation

This host has no AArch64 cross-assembler or compiler, QEMU AArch64,
`armasm64`, LLVM ARM assembler, Morello/CHERI SDK, CHERI emulator, or CHERI
hardware. Consequently this phase does not claim ARM64 assembly/execution or
CHERI assembly/execution evidence.

When those targets are available, the remaining gate is:

1. cross-assemble the generated ELF/COFF/Mach-O/Win64 outputs with their
   supported AArch64 toolchains;
2. assemble the CHERI output with a Morello-capable toolchain;
3. run `InterfaceTests.CtInverseMod256` and
   `InterfaceTests.CtxInverseMod384` on native AArch64 and Morello or the
   instruction emulator; and
4. confirm the deterministic boundary inputs complete without capability
   faults and satisfy the multiplicative identity.

This toolchain gap is a release-evidence requirement, not an open defect in
the reviewed generator change.
