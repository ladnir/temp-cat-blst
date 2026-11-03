# File Review: assembly integration and generated artifacts

## Record

| Field | Value |
| --- | --- |
| Paths | `blst/generated/assembly.S`, `scripts/build/refresh.sh`, `scripts/build/generate-win64-exports.sh`, `blst/generated/win64/blst.def`, and supported checked-in assembly outputs |
| Roadmap phase / queue | Phase 13 assembly trust-chain integration closure |
| Reviewed commit | `5fbe3a5a732264727087f8dcc4db172d80f87f14` |
| Reviewer | Codex audit team |
| Date | 2026-08-13 |
| Disposition | `generated_production` and supporting generators |
| Subsystem | Assembly selection, generated-object integration, and Windows export metadata |
| Profiles | x86-64 ELF/COFF/Mach-O/MASM; supported Apple ARM64 Mach-O; static packaging metadata |
| Generation/provenance | Audited BLST-derived generators and translators with Category-blst deterministic orchestration |
| Upstream references | Generator-specific references are recorded in `docs/audit/ASSEMBLY_AUDIT_QUEUE.md` |
| Related files | `blst/CMakeLists.txt`, `blst/common/cpuid.cpp`, `bindings/blst.h`, `blst/generated/assembly-output-manifest.txt`, `scripts/build/verify-assembly.sh`, `scripts/audit/check_assembly_translator_objects.sh` |

## Purpose and trust boundary

`assembly.S` selects checked-in providers from public target and capability
macros. It contains no arithmetic. The selected generated files implement the
raw arithmetic and SHA contracts reviewed in Batches B through E.

The export generator preprocesses the supported C umbrella header and emits
the future Windows DLL export allowlist. Shared-library construction remains
explicitly disabled, so the `.def` file is packaging metadata rather than an
active Release linker input.

`refresh.sh` is a bounded operation dispatcher for assembly, Windows exports,
and Rust bindings. It accepts an explicit operation and exact argument count,
then replaces itself with the selected generator.

## Contracts and invariants

- Assembly is available only on supported 64-bit targets. An unsupported
  assembly target fails preprocessing. A build that disables assembly may
  intentionally preprocess an empty aggregate translation unit.
- x86 selection includes the target object format and the provider families
  required by `BLST_ADX` and `BLST_PORTABLE`. Public capability dispatch occurs
  inside the selected owners.
- Apple ARM64 selects only the audited Mach-O ARM providers. Non-Apple ARM64
  and CHERI outputs remain retained experimental artifacts.
- Every checked-in artifact must be reproduced byte-for-byte by its audited
  generator and pinned translator. The allowlist must match the checked tree.
- The Windows export generator resolves repository-root includes, writes
  through a private temporary file, and atomically replaces only its explicit
  output path.
- `blst.def` contains each public umbrella-header symbol exactly once. Every
  listed symbol must exist in the corresponding complete static archive before
  future shared-library work consumes the file.
- Platform and capability controls are public. This layer handles no secret
  values and adds no secret-dependent flow.

## Dependency and call-path map

Non-MSVC CMake builds compile `assembly.S`. MSVC builds select the corresponding
checked MASM files directly. `assembly.S` includes generated ELF, COFF, or
Mach-O text according to compiler target macros. The resulting symbols satisfy
the raw declarations in `backend_asm.h` and are selected by `backend.h`.

`generate-assembly.sh` owns deterministic artifact production;
`verify-assembly.sh` regenerates twice and compares the complete allowlist.
`check_assembly_translator_objects.sh` assembles aggregate target objects and
compares their defined external symbols with source declarations. The public
C umbrella header supplies the input to `generate-win64-exports.sh`.

The deleted zero-byte `build/assembly.h` had no caller. The deleted
`blst/generated/win64/dll.c` belonged to the removed batch-file DLL path. Its
initialization and conforming memory helpers already live in `cpuid.cpp` behind
`BLST_DLL_MAIN` for a future audited shared-library implementation.

## Upstream/specification comparison

Arithmetic and translator deltas are recorded per owner in the assembly queue.
This packet changes only integration policy. Category-blst rejects 32-bit
targets, uses split C headers under the compatibility umbrella, and currently
rejects shared builds. The aggregate selector and export metadata now match
those local policies.

## Checklist evidence

### 1. Identity and purpose

The aggregate selector has Apache-2.0, BLST, and Category Labs attribution.
Generator, selector, export, and dispatcher ownership are distinct.

### 2. Contracts and invariants

Target, capability, artifact, export, failure, timing, and ownership contracts
are stated above. Raw arithmetic contracts remain at `backend_asm.h`.

### 3. Upstream and specification lineage

Every arithmetic owner is bound to its upstream row and exact generated corpus
through `ASSEMBLY_AUDIT_QUEUE.md`. Local orchestration changes are separately
recorded as `AUD-359` through `AUD-378` and `AUD-380` through `AUD-382`.

### 4. Correctness

The arithmetic owners passed their independent byte, field, SHA, and alias
oracles. Exact regeneration and target-object symbol comparisons bind those
results to the checked outputs. This integration packet changes no generated
arithmetic byte.

### 5. Memory and lifetime safety

The selector owns no runtime storage. Export generation uses a guarded private
temporary directory and an explicit output. The obsolete DLL memory helpers
with nonstandard parameter types were deleted.

### 6. Undefined and implementation-defined behavior

Target selection relies only on compiler platform macros. Unsupported 32-bit
assembly no longer silently produces an empty object. Object-format, unwind,
and ABI checks cover the supported generated dialects.

### 7. Side channels and cryptographic hygiene

Selection depends only on public build and CPU capability values. Arithmetic
constant-flow evidence is owned by the source generators and code-generation
reviews. This packet adds no runtime instruction.

### 8. Randomness and secret lifecycle

N/A: generation and selection are deterministic and process no secrets.

### 9. API and ABI

The checked export file now contains all 342 public C names, no duplicate, and
no removed constant. Every name appears in the current x86-64 Release archive.
Shared-library ABI construction remains intentionally disabled.

### 10. Concurrency

Runtime code is immutable. Generators use independent temporary directories;
the checked-tree verifier does not mutate tracked outputs.

### 11. Resource exhaustion and denial of service

Generation processes a fixed repository corpus. CI jobs have a 15-minute
timeout. Runtime assembly schedules and stack frames are reviewed per owner.

### 12. Portability and build behavior

Aggregate preprocessing passes for supported x86-64 and deliberate no-assembly
profiles. An enabled unsupported 32-bit profile fails. Apple ARM64 has prior
hosted execution. Experimental targets retain their separate promotion gates.

### 13. Performance-sensitive structure

No provider body, dispatch instruction, call boundary, or generated byte
changed. The export and CI changes are build-time only.

### 14. Tests and evidence

The focused checks and existing batch evidence are summarized below. No broad
build or benchmark was required.

### Specialized additions

The generated artifact review follows the owning source and translator review.
A generated file must be reopened whenever its owner, translator, output
format, or exact bytes change.

## Findings

`AUD-009`, `AUD-041`, `AUD-359` through `AUD-378`, and `AUD-380` through
`AUD-382` as applicable to each owner and artifact family.

## Validation

- All 97 retained assembly outputs reproduce deterministically.
- Aggregate ELF, COFF, Mach-O, MASM, and supported Apple ARM64 object symbol
  sets match their source declarations; applicable unwind records decode.
- Independent Batches B-E arithmetic and SHA evidence remains attached to the
  unchanged generated bytes.
- Supported x86-64 aggregate preprocessing: passed.
- Assembly-disabled unknown-target preprocessing: passed with empty output.
- Enabled 32-bit target preprocessing: rejected as required.
- Windows export generation: 342 unique names; exact regeneration passed.
- Existing x86-64 Release archive: all 342 export names discovered; none
  missing.
- Workflow and generator shell syntax: passed.
- No shared-library link, product rebuild, benchmark, sanitizer, CUDA,
  portable-AVX Release, BTx/BTE audit, or remote job ran.

## Automated evidence

| Tool/profile | Result | Boundary |
| --- | --- | --- |
| Deterministic generator verification | 97/97 exact outputs | Generator and translator execution |
| Aggregate object inspection | Exact public symbols and decoded supported unwind | Does not execute target arithmetic |
| Independent arithmetic tests | Passed per Batches B-E | Platform execution gaps remain per finding |
| Audit-integrity export check | Added; local equivalent passed | Hosted execution awaits the pushed revision |
| Sanitizers and fuzzing | Not rerun | Milestone evidence remains separate |

## Residual risk and decision

The supported generated owners and integration paths are complete for
file-by-file review. Shared-library construction remains disabled and requires
its own packaging and ABI packet. The export list is ready evidence, not a
claim that a DLL was linked. Experimental ARM64 and CHERI artifacts remain
outside the supported completion set. Native platform, sanitizer, fuzz, and
release-matrix gaps remain recorded at their owning findings.
