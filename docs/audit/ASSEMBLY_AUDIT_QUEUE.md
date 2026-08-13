# Scalar Assembly Trust-Chain Audit Queue

## Freeze and scope

The assembly audit resumes from commit
`aa7c2a777dd82cbc4ac5254bbb996e077ac98939`. It covers the 25 handwritten
Perl owners under `blst/asm`, the generation/verification control files under
`build`, and every checked-in generated output named by
`build/assembly-output-manifest.txt`.

The manifest names 100 outputs: 23 COFF, 23 ELF, 23 Mach-O, 22 Win64, and 9
CHERI files. Generated files are reviewed with their owning generator and
translator; deterministic byte reproduction alone does not establish
arithmetic, ABI, constant-time, or memory-safety correctness.

Benchmarks remain permission-gated and sequential. This queue does not build
CUDA, audit BTx/BTE, run a broad product suite, or treat portable AVX as native
deployment evidence.

## Dependency order

| Batch | Handwritten owners | Generated outputs | Primary questions |
| --- | --- | ---: | --- |
| A. Translation and generation boundary | `arm-xlate.pl`, `x86_64-xlate.pl`, `generate-assembly.sh`, `verify-assembly.sh`, output allowlist | all 100 structurally | provenance, supported flavours, deterministic environment, error propagation, output completeness, ABI/unwind translation |
| B. Add/subtract and division helpers | five `add_mod*` owners and two `div3w*` owners | 31 | exact modular ranges, borrow/carry behavior, aliasing, count bounds, constant-time control flow |
| C. Montgomery multiplication | two ARM `mul_mont*` and four x86 q/x `mul*_mont*` owners | 26 | reduction bounds, ADX versus generic equivalence, Fp/Fp2 layout, repeated-square loops, aliasing |
| D. Inversion and square testing | seven `ct*inverse*`/`ct_is_square*` owners | 31 | fixed iteration schedules, scratch bounds, zero behavior, sign/correction steps, CHERI capability bounds |
| E. SHA-256 | ARM, x86 SHA-extension, and portable-x86 owners | 12 | block-count preconditions, dispatch symbols, endian loads, SHA-extension fallback, constant-time behavior |

The two translators are reviewed before arithmetic owners because every output
passes through one of them. Arithmetic batches may proceed independently only
after Batch A fixes and evidence are integrated.

## Batch A inventory

- `x86_64-xlate.pl` is byte-identical to pinned upstream BLST commit
  `54e6e55674722fc2797ebb4bbb71b26d881eb4b8`, apart from the repository path
  and executable mode.
- `arm-xlate.pl` differs from that snapshot only in the reviewed CHERI
  `scbnds` support and multi-digit capability-register translations introduced
  by `ce8e666`. The focused CHERI operand probe passes at the freeze commit.
- The generation wrapper normalizes CRLF input and output, fixes locale and
  timezone, requires an empty destination, and checks the exact path allowlist.
- The verification wrapper generates twice, compares both generations, then
  compares every normalized result with the checked-in tree.
- Shell syntax checks pass for both wrappers. One focused ARM generator emits
  translated Linux assembly successfully through the stdout convention.

## Batch A dispositions

1. `AUD-359` retains all outputs but classifies the 36 non-Apple ARM64 and
   CHERI artifacts as experimental. The 9 Apple ARM64 and 55 x86-64 outputs
   remain in the production claim.
2. `AUD-361` makes all 23 arithmetic owners propagate translator completion
   failure and makes the wrapper reject empty raw output. Focused failure
   injection covers one x86 and one ARM owner; a structural check covers the
   complete owner set.
3. `AUD-362` preserves the pinned generic x86 translator but invokes it from
   checked generation with `CC=false`. This deterministically selects the
   empty object-format prefix used by the checked-in x86-64 artifacts.
4. `AUD-360` changes the experimental MSVC ARM64 source selection from pointer
   width to the root configuration's validated target architecture. The local
   MSVC installation lacks the ARM64 toolset, so this finding remains fixed
   rather than validated by a native target configuration.

The generation-boundary changes are committed at `61d3b66`. They alter no
generated output bytes and do not promote any experimental target.

## Batch A validation plan

- compare both translators against the pinned upstream snapshot and record
  every local divergence;
- exercise a compact directive/operand corpus for every generated flavour;
- inject translator and output failures and require nonzero generator status;
- run the existing CHERI operand probe and Win64 symbol allowlist;
- regenerate twice and compare all 100 outputs with the checked-in tree;
- inspect representative Mach-O, ELF, COFF, and Win64 symbols/unwind metadata;
- compile only small representative output objects where a local assembler is
  available; and
- run no benchmark or broad product test.

## Batch A closure

Batch A is complete at `1642ff4` for Category-blst's checked generator corpus.
This is not a claim that the inherited translators correctly accept arbitrary
OpenSSL perlasm input. Perl syntax checks, the CHERI operand probe, translator
and output-failure injection, paths containing spaces, and exact
two-generation verification all pass. Verification reproduces every one of
the 100 checked outputs byte-for-byte. A supported MSVC x86-64 Release
configuration still selects the x86 family.

The complete aggregate source cross-assembles with Clang 18 as x86-64 ELF,
x86-64 COFF, x86-64 Mach-O, Apple ARM64 Mach-O, experimental ARM64 ELF, and
experimental ARM64 COFF. For every object, the defined external symbol set
exactly matches the preprocessed source declarations. ELF exposes decoded
`.eh_frame`; COFF exposes decoded `.pdata` and `.xdata`; Mach-O applies the
expected underscore decoration. Existing MSVC/MASM objects retain the exact
38-entry portable fallback surface, and the representative SHA object has
decoded SEH records for its saved general and XMM registers.

Native ARMASM acceptance, CHERI assembly/execution, and non-Apple ARM runtime
qualification remain target-promotion evidence rather than Batch A gaps.
Arithmetic correctness, ABI behavior inside each routine, constant-time
control flow, and execution proceed in Batches B through E.

## Batch B inventory and findings

Batch B covers five add/subtract owners, two division owners, and 31 generated
outputs. The public field layers call the selected wrappers with reduced Fr,
Fp, or Fp2 values. Fp12 uses the double-width routines before Montgomery
reduction. Scalar and AVX group multiplication use the division helpers to
split public or canonical scalars by the BLS parameter.

The arithmetic inspection found no confirmed wrong result. It opened four
findings that prevent closure:

1. `AUD-365` records the zero-length prefetch range underflow.
2. `AUD-366` records 32-bit count truncation in four raw x86 shift loops.
3. `AUD-367` requires independent boundary and aliasing oracles.
4. `AUD-368` requires exact range, modulus, and division-estimate contracts.

The narrow arithmetic contracts are intentional. Single-width add, subtract,
shift, and small-multiplier routines accept reduced operands. Double-width
add and subtract accept operands below `p * 2^384` and reduce modulo that
boundary. `div_3_limbs` accepts a normalized divisor and a dividend prefix
strictly below that divisor. It returns an exact or one-high quotient estimate.
The 128-bit remainder helper corrects one-high estimates. The 64-bit helper
requires an exact estimate.

Twenty focused x86-64 assembly tests passed before remediation. A separate
integer model checked 200,000 normalized strict-bound division estimates.
No build, benchmark, sanitizer sweep, CUDA test, or BTx/BTE work ran during
the read-only inventory.

## Batch B remediation evidence

Commit `2333650` implements the approved boundary repairs. The selected
prefetch function returns for a zero length. Positive fixed sizes still compile
to a direct tail call to `vec_prefetch_asm`. Variable sizes add one branch on
the public length. The raw assembly contract remains positive-length only.

All four raw x86 shift loops now decrement the complete 64-bit `size_t` count.
The loop bodies and fixed-count selected code remain unchanged. Deterministic
generation reproduces all 100 checked outputs after the intentional 13-file
delta. Structural inspection finds no 32-bit decrement in the affected owners
or outputs.

The new byte-level oracle does not call portable limb arithmetic. It covers:

- 256-bit and 384-bit modular addition, subtraction, and shifts;
- zero, one, `p-1`, and deterministic dense reduced operands;
- exact output aliasing with either binary input;
- addition and subtraction modulo `p * 2^384`;
- canonical scalar checks and checked add/subtract boundaries; and
- quotient estimation, one-high correction, saturation, and exact 64-bit
  division.

Focused MSVC 19.50 x86-64 Release execution passes 6/6. A Debug compile of the
widened Fp12 path accepts the new operand and alias assertions. The optimized
selector probe compiles under `/O2 /Ob3 /DNDEBUG /W4 /WX`; its object hash is
`C99F5717D45703F738134FA65C73FEB61A3BF57B0FCEA63A7C3529E2D3819C58`.

`AUD-366` and `AUD-368` are validated. `AUD-365` and `AUD-367` are fixed but
retain Apple ARM64 execution as supported-platform evidence. No benchmark,
sanitizer, CUDA, portable-AVX Release, full test sweep, or BTx/BTE work ran.
