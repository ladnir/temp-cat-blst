# Scalar Assembly Trust-Chain Audit Queue

## Freeze and scope

The assembly audit resumes from commit
`aa7c2a777dd82cbc4ac5254bbb996e077ac98939`. It began with 25 handwritten
Perl owners under `blst/asm`, the generation/verification control files under
`build`, and every checked-in generated output named by
`build/assembly-output-manifest.txt`.

Batch E removed one unreachable scalar x86 SHA owner. The retained tree has 24
handwritten Perl owners.

The manifest now names 97 outputs: 22 COFF, 22 ELF, 22 Mach-O, 22 Win64, and
9 CHERI files. Generated files are reviewed with their owning generator and
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
| E. SHA-256 | ARM and x86 SHA owners | 9 | block-count preconditions, dispatch symbols, endian loads, hardware selection, no-ISA fallback, constant-time behavior |

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

`AUD-365` through `AUD-368` are validated. Apple scalar-assembly run
`31712201574` supplies the supported Apple evidence for the shared owners. No benchmark,
sanitizer, CUDA, portable-AVX Release, full test sweep, or BTx/BTE work ran.

## Batch C inventory and findings

Batch C covers six Montgomery multiplication owners and 26 generated outputs.
The generators remain arithmetic-equivalent to pinned upstream. Local deltas
are limited to the approved translator boundary, capability-symbol spelling,
and full-width public repeated-square counters. No wrong-result defect was
confirmed during source inspection.

The raw 256-bit sparse multiplier has an asymmetric domain. Its first operand
is reduced modulo r, while its second operand may occupy all 256 bits. Sparse
squaring requires a reduced operand. REDC accepts T only when `T < rR` or
`T < pR`. The 383-bit repeated-square routine requires a positive public count,
permits output reuse of its first operand, and forbids output reuse of its
final multiplier. Fp2 assembly assumes contiguous real and imaginary Fp
components. The 382x routines return two contiguous double-width components.

The review opened three findings:

1. `AUD-369` records correlated tests and two invalid test-input domains.
2. `AUD-370` records missing raw Montgomery contracts.
3. `AUD-371` records unpinned C++ layouts at the assembly boundary.

The approved repair uses a branchy byte-level model. For Montgomery output y,
the model checks the unique reduced solution through
`yR mod m = ab mod m`. REDC and from-Montgomery use the corresponding
congruence. Fp2 and 382x checks derive both component equations independently.
No production arithmetic, generator instruction, runtime check, or benchmark
change is part of this packet.

## Batch C remediation evidence

Commit `a29895f` replaces closure comparisons with an independent branchy
byte model. The oracle uses schoolbook byte multiplication and bitwise modular
division. It does not call portable limb or Montgomery arithmetic. It covers:

- selected, forced-generic, and direct ADX 256-bit and 384-bit providers;
- the sparse multiplier's full-width second operand;
- reduced sparse squares and valid `T < mR` REDC inputs;
- from-Montgomery over the complete radix width;
- exact output aliasing at every supported boundary;
- Fp2 real and imaginary component equations;
- unreduced 382x component equations; and
- repeated-square counts 0, 1, 2, and 5, including the optimized 383-bit path.

The generic x86 test clears only the public ADX capability bit for its scope,
then restores the process capability word. This forces the q implementation
on ADX hardware. The direct x implementation remains separately exercised.

Focused MSVC 19.50 x86-64 Release execution passes 15/15 assembly-audit tests.
Eleven retained Montgomery interface tests also pass after their invalid input
generators were corrected. Compile-time assertions pin all six assembly-visible
field shapes. No generated assembly or production arithmetic changed.

`AUD-370` is validated. `AUD-369` and `AUD-371` are fixed pending supported
Apple ARM64 execution and compilation. No full build, benchmark, sanitizer,
CUDA, portable-AVX Release, or BTx/BTE work ran.

## Batch D inventory and findings

Batch D covers seven fixed-schedule inverse and quadratic-residue owners and
31 generated outputs. All seven generators retain upstream arithmetic. Local
differences before this packet were the approved fail-closed translator calls
and the public x86 capability-symbol spelling.

The schedules are fixed:

- the 256-bit inverse performs `15 * 31 + 47 = 512` divsteps;
- the generic 384-bit inverse performs `12 * 62 + 24 = 768` divsteps;
- the ADX 384-bit inverse performs `23 * 31 + 55 = 768` divsteps; and
- the residue test performs `24 * 30 + 48 = 768` divsteps.

Secret-dependent choices use masks or conditional moves. Branches depend only
on fixed loop counters. The additional x86 branch selects a provider from the
public capability word. Each inverse uses a fixed stack frame with one aligned
512-byte scratch region. The residue test uses one aligned 256-byte region.

The review opened three findings. `AUD-372` records inaccurate raw contracts.
`AUD-373` records correlated correctness evidence. `AUD-374` confines the
experimental CHERI residue-test scratch capability to its intended region.
No arithmetic defect was confirmed.

## Batch D remediation evidence

Commit `83c11ba` states the exact raw contract. For width `w` and `R=2^w`, a
nonzero invertible input produces a bounded double-width value `T` satisfying
`T * inp = R^2 (mod mod)`. Zero produces zero. Field callers normalize this
private representation. The residue test accepts a reduced value modulo an
odd prime and returns true for zero and nonzero quadratic residues.

The independent oracle uses branchy byte multiplication, bitwise modular
division, and binary modular exponentiation. It does not call portable inverse
or residue-test code. Focused MSVC 19.50 x86-64 Release execution passes all
five new tests for selected, portable, direct ADX, and forced-generic inverse
providers plus assembly and portable residue tests. Three retained interface
tests also pass.

The CHERI ARM generator now bounds its scratch capability to 256 bytes and
derives each XOR-selected peer from a bounded capability. Deterministic
verification reproduces all 100 checked outputs. The CHERI translator fixture
passes. Preprocessing the supported Apple ARM64 path before and after the
change produces identical ordinary instructions.

`AUD-372` is validated. `AUD-373` remains fixed pending supported Apple
execution of the independent oracle. `AUD-374` remains fixed because CHERI is
experimental and no Morello toolchain is available. No benchmark, broad test
sweep, sanitizer, CUDA, portable-AVX Release, or source change under BTx/BTE
ran. Batch E, SHA-256, is next.

## Batch E inventory and findings

Batch E covers the ARM and x86 SHA-256 owners and nine generated outputs. The
compression arithmetic retains the upstream SHA-extension, SSSE3, and NEON
implementations. The inspection confirmed one dispatch defect: x86 selected
SSSE3 whenever SHA-NI was unavailable, even when the CPU lacked SSSE3. The
separate scalar x86 generator was unreachable because its only guard,
`BLST_PORTABLE_2`, was never defined.

Four findings define the remediation:

1. `AUD-375` requires a safe x86 provider chain for every supported CPU.
2. `AUD-376` confines the experimental CHERI message schedule capability.
3. `AUD-377` states the raw SHA block, alias, and timing contracts.
4. `AUD-378` requires provider-specific known-answer evidence.

## Batch E remediation evidence

Commit `2f59883` gives x86 SHA compression three explicit providers. SHA-NI
wins when CPUID leaf 7 advertises it. SSSE3 is selected only when CPUID leaf 1
advertises SSSE3. The portable C compression routine is the no-ISA fallback.
`BLST_PORTABLE` continues to govern generic arithmetic assembly; it no longer
implies a SHA fallback policy. SGX-LVI builds preserve the inherited policy of
not selecting SHA-NI.

The public dispatcher is stackless. Its normal SHA-NI path performs one public
capability test and a direct tail jump. SSSE3 adds one test. The fallback adds
no call frame, allocation, conversion, or message-dependent access. The dead
scalar x86 generator and its three generated outputs are removed, reducing the
checked corpus from 100 to 97 files.

The ARM owner retains its existing public capability dispatch: hardware SHA
instructions are preferred, and the generic NEON implementation remains the
fallback. CHERI schedule pointers are now bounded to the 64-byte schedule
region. Preprocessing removes those bounds from ordinary Apple ARM64 code. The
nonblank Apple instruction projection is byte-identical before and after the
change, with SHA-256 hash
`93565BB099430E735001EB1BF968B4D0AC51C617FC1B90912C08B345716EAA0D`.

Forced x86 C, SSSE3, and available SHA-NI providers pass the `abc` and
two-block deterministic known-answer vectors. Six focused SHA tests pass in
the supported MSVC Release configuration. Five pass in the nonportable ADX
configuration. Deterministic generation reproduces all 97 outputs, and the
aggregate x86-64 and ARM64 translator-object checks pass. The prior supported
Apple Silicon run remains applicable because the Apple instruction projection
did not change.

`AUD-375`, `AUD-377`, and `AUD-378` are validated. `AUD-376` is fixed; CHERI
execution remains an experimental platform-promotion requirement. No
benchmark, broad test sweep, sanitizer, CUDA, portable-AVX Release, or BTx/BTE
work ran. Batch E completes the handwritten scalar assembly-family audit.

## Batch F integration closure

Commit `5fbe3a5a` closes the supported aggregate and generated-artifact layer.
The aggregate selector now rejects unsupported enabled targets, while an
explicit no-assembly build may still emit an empty translation unit. The
Windows export generator resolves repository headers from a clean checkout and
reproduces all 342 public C symbols. Every symbol exists in the current Release
archive, and audit-integrity CI now rejects drift in the checked export list.

The 52 retained x86-64 outputs and nine Apple ARM64 Mach-O outputs inherit the exact
generation, object, arithmetic, and platform evidence from Batches A through
E. They are complete for file review. The 36 non-Apple ARM64 and CHERI outputs
remain experimental and retain their platform-promotion gates. The retained
corpus remains 97 files because this packet changes no generated byte.

The obsolete Windows DLL owner and zero-byte assembly header were deleted.
Shared-library construction remains disabled and requires a separate ABI
packet. No product rebuild, benchmark, sanitizer, CUDA, portable-AVX Release,
BTx/BTE audit, or remote job ran.
