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

## Questions requiring disposition

1. The manifest currently classifies all 100 outputs as production. Thirty-six
   are non-Apple ARM64 or CHERI artifacts, although those platforms are
   experimental. The likely repair is a manifest-only scope correction:
   retain the generators and artifacts but make only 9 Apple ARM64 and 55
   x86-64 outputs part of the production claim.
2. The 23 arithmetic generators generally open translator pipes and close
   `STDOUT` without uniformly checking failure. `verify-assembly.sh` catches a
   mismatched checked-in tree, but `generate-assembly.sh` can provide a weaker
   success boundary. The likely repair is mechanical error propagation with
   zero generated-byte changes.
3. The upstream x86 translator probes `CC` through a shell command for the
   MinGW symbol prefix. Category-blst's checked x86-64 profiles all use the
   empty prefix, but generation can still depend on the caller environment.
   Decide whether the wrapper should force the known deterministic prefix or
   preserve the generic upstream probe and document it.

These are audit candidates, not approved product changes. Batch A will assign
finding IDs only after their impact and smallest fixes are agreed.

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

