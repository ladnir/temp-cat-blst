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

## Batch A evidence and remaining work

At `61d3b66`, Perl syntax checks, the CHERI operand probe, translator-failure
injection, and exact two-generation verification all pass. Verification
reproduces every one of the 100 checked outputs byte-for-byte. A supported
MSVC x86-64 Release configuration still selects the x86 family.

The two translators remain in review for the complete directive/flavour corpus
and representative object symbol/unwind inspection. Those tasks, plus native
target qualification, are intentionally distinct from the now-validated
generation failure boundary. Arithmetic correctness, ABI, constant-time, and
target execution remain assigned to Batches B through E.
