# Phase 3A Assembly Reproducibility Evidence

## Record

- Date: 2026-07-30
- Baseline commit: `6103816`
- Implementation commit: `942404d`
- State: validated and committed
- Scope: `AUD-009` and review-discovered `AUD-016`
- Benchmarks: disabled and not run
- CUDA: disabled

This phase establishes provenance for the existing fork assembly. It does not
import any upstream instruction-schedule, symbol-visibility, or translator
change. Those remain Phase 3B/3C work.

## Generator inventory

The checked-in generator set contains:

- 14 x86-64 code generators, plus `x86_64-xlate.pl`;
- 9 ARMv8 code generators, plus `arm-xlate.pl`;
- 55 x86-64 outputs: ELF, COFF, Mach-O, and 13 MASM files because the portable
  SHA-256 generator intentionally has no MASM output;
- 45 ARMv8 outputs: ELF, COFF, Mach-O, MASM, and CHERI;
- 100 generated assembly files in total.

`blst/generated/assembly-output-manifest.txt` is the exact output allowlist. Generation
refuses a nonempty output directory and never writes directly into the
checked-in output tree.

## Reproducibility repair

`scripts/build/refresh.sh` previously searched the removed `src/asm` directory and
coupled assembly, Win64 export-definition, and Rust-bindgen refreshes. It now
dispatches those operations independently:

- `generate-assembly`;
- `verify-assembly`;
- `generate-win64-exports`;
- `generate-rust-bindings`.

Assembly generation copies the checked-in Perl sources into a temporary
canonical-LF source tree before execution and normalizes generated text to LF.
This prevents Windows `core.autocrlf` from changing embedded here-document
bytes and therefore changing generated output.

The checked-in Win64 outputs for `mulq_mont_256` and `mulq_mont_384` contained
manual `BLST_PORTABLE` substitutions while their generators emitted
`__BLST_PORTABLE__`. The generators now encode that existing MASM-only fork
spelling explicitly. Phase 3B will replace this temporary compatibility rule
with the final unified MASM macro policy.

Phase 3B subsequently removed this temporary rule and validated the unified
policy; see `PHASE_03B.md`.

Repository attributes require LF for generator scripts, shell scripts, the
allowlist, and checked-in assembly. Verification still normalizes an existing
Windows checkout before comparison so it works before a fresh checkout.

## Validation

### Generator checks

`scripts/build/refresh.sh verify-assembly` generated two independent temporary trees.
Every corresponding output was byte-identical, both trees matched the
100-entry allowlist, and all outputs matched the checked-in tree after
repository newline normalization.

The split Win64 export generator produced byte-identical output twice. The
Rust-bindgen generator was unavailable on the validation host; its preflight
was tested with a missing executable and failed without creating an output
file. Assembly verification does not depend on bindgen.

Shell syntax validation passed for every new or changed generator wrapper.
`git diff --check` passed, and no checked-in generated assembly changed.

### Build and test matrix

| Profile | Configuration | Result |
| --- | --- | --- |
| Scalar assembly | MSVC Debug, default scalar API, MASM/ADX enabled, CUDA/benchmarks off | Build passed; 310/310 registered tests passed |
| AVX-512 portable | MSVC Debug, default AVX API, portable AVX-512 plus MASM/ADX, CUDA/benchmarks off | Build passed; 142/142 registered tests passed |

The AVX-512 portable run included the retained FFT and middle-product
algorithm tests. No benchmark executable was built or run.

## Review-discovered AUD-016

The newly isolated Win64 export generator is deterministic, but its output
does not match `blst/generated/win64/blst.def`. The checked-in file contains 233 export
entries, while the current supported C header generates 318. No tracked build
currently consumes the definition file.

This is recorded as `AUD-016` rather than silently replacing the Windows ABI
surface during assembly work. It should be resolved with the supported C ABI,
shared-library packaging, CatBLST naming, and configurable prefix policy.

Later resolution: the complete current 342-symbol compatibility surface now
matches a fresh deterministic generator run byte-for-byte. Shared packaging
and configurable prefixing remain separate Phase 20 work.

## Residual validation

- ARMv8 and CHERI outputs match their checked-in baselines but were not
  assembled or executed on this Windows x86-64 host.
- Cross-version Perl and non-WSL host reproduction should be added to CI.
- Rust bindings were not regenerated because bindgen is not installed.
- `AUD-010` MASM symbol visibility and the final portable macro policy were
  subsequently validated in Phase 3B.
- Accepted upstream generator and translator changes remain Phase 3C work and
  must be applied only through this verified path.
