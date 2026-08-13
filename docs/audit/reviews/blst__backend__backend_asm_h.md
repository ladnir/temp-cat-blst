# File Review: `blst/backend/backend_asm.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/backend/backend_asm.h` |
| Roadmap phase / queue | Phase 13, scalar assembly interface |
| Reviewed commit | `2f59883caf67191aae95abb905598f405f8139e8` |
| Reviewer | Codex audit team |
| Date | 2026-08-13 |
| Disposition | `in_scope_production` |
| Subsystem | Raw scalar assembly declarations and provider availability |
| Profiles | x86-64 ADX, x86-64 q, combined portable x86-64, Apple/AArch64 scalar |
| Generation/provenance | Locally composed declaration header with upstream BLST `src/vect.h` lineage |
| Upstream references | Upstream snapshot `54e6e55674722fc2797ebb4bbb71b26d881eb4b8`; local provider repair `33914f2`; generated assembly tree |
| Related files | `backend.h`, `common/vect.h`, `build/assembly.S`, generated x86-64/AArch64 outputs, `test/Interface_tests.cpp` |

## Purpose and trust boundary

This internal header is the C-linkage contract between C++ scalar code and
handwritten/generated assembly. It declares fixed-width modular arithmetic,
wide Fp2 operations, inverse and division kernels, vector predicates, and
SHA-256 helpers. The declarations do not select a production backend; selection
belongs to the forced-inline wrappers in `backend.h` and the scalar field layer.

Operands can contain secret scalars, field elements, and key-derivation state.
The interface therefore requires exact ABI widths, correct const direction,
no allocation or exception behavior, and an architecture/provider surface that
matches what the build actually links.

## Contracts and invariants

- CatBLST uses 64-bit limbs exclusively.
- The raw `asm_*` annotations describe byte/limb shapes while intentionally
  decaying to pointers at the assembly ABI. They create no C++ typed access and
  therefore no strict-aliasing mismatch.
- Raw vector predicates require a positive byte count divisible by 16.
- Raw SHA compression and byte-copy functions require positive counts.
  Forced-inline C++ boundaries define zero-count behavior before assembly.
- Single-width modular operations require reduced operands. Exact output
  aliasing with either input is supported; partial overlap is not.
- Double-width add/subtract operates modulo `p * 2^384` and requires each
  operand below that boundary.
- Division routines require their documented normalized divisor and exact or
  one-high quotient estimate.
- x86-64 portable assembly builds provide both x and q public providers.
- Non-portable GNU/Clang x86-64 builds provide the configured x or q provider.
- Win64 x86-64 builds compile both public MASM providers.
- AArch64 builds provide only generic provider names.
- Provider names are declarations, never preprocessor dispatch aliases.
- The header owns no state and performs no runtime CPU dispatch.

## Provider matrix

| Target/configuration | x declarations | q declarations |
| --- | --- | --- |
| x86-64, `BLST_PORTABLE=ON` | yes | yes |
| x86-64 GNU/Clang, ADX non-portable | yes | no |
| x86-64 GNU/Clang, q non-portable | no | yes |
| MSVC x64 | yes | yes |
| Apple/AArch64 | no | yes |

The 16 x-provider and 16 q-provider declarations cover 256-bit Montgomery
operations, 384-bit/Fp2 Montgomery operations, inverse, sign/parity, repeated
squaring, and 382-bit wide multiplication/squaring. Forty-three common
declarations cover fixed modular operations, vector helpers, SHA, division,
and wide add/subtract. The optional SHA-extension declaration and internal C
fallback bring the exact inventory to 77.

## Checklist evidence

### Identity, history, and provenance

The file, its complete history, all 77 declarations, generated providers,
direct production/test callers, and upstream declaration shapes were inspected.
The header now carries Apache-2.0/SPDX and explicit local/upstream provenance.

### Correctness and ABI

Every declaration was matched to a generated symbol and its documented operand
shape. In particular, 256- and 384-bit constant-time inverse outputs are the
expected double-width 512/768-bit values, and 382-bit Fp2 outputs are two
768-bit values. Input constness and scalar return widths match the generated
contracts. The C symbol names and calling conventions are unchanged.

### Memory safety and undefined behavior

The header performs no accesses itself. `AUD-056` records and fixes the former
failure to enforce positive raw counts at selected vector/SHA boundaries. The
`void*` raw ABI annotations do not dereference through an incompatible C++
type. Fixed shape enforcement remains the responsibility of audited callers and
their concrete field representations.

### Side channels and performance

No data-dependent selection was introduced. Provider availability is resolved
by preprocessing; production wrappers call raw symbols directly. MSVC `/O2
/Ob3` fixed-count probes tail-call `vec_is_zero_16x`,
`blst_sha256_bcopy`, and `blst_sha256_block_data_order` with no zero branch,
copy, conversion, or intermediate call.

### Portability and build behavior

`AUD-057` removes the persistent generic-to-ADX and AArch64 x-to-generic macro
aliases. Header order can no longer rename subsequent code. The matrix mirrors
`build/assembly.S` and the Win64 CMake object selection. Local MSVC combined
and non-portable configurations compile and retain their exact fallback symbol
allowlists. Supported run `31665323681` builds and tests the direct generic
provider on Apple Silicon.

### API, exceptions, concurrency, and resources

This is an internal native C++/assembly boundary, not a public C ABI. It
allocates nothing, throws nothing, stores no pointers, owns no mutable global or
thread-local state, and creates no worker or synchronization behavior.

## Findings

- `AUD-041` — validated: scalar provider selection/build configurations.
- `AUD-056` — validated: zero-length selected boundaries no longer enter raw
  positive-count assembly loops.
- `AUD-057` — validated: provider aliases are removed, direct differential
  evidence is repaired, and supported Apple/AArch64 execution passes.

## Focused validation

- Combined Win64 portable configuration compiled the exact changed
  `Interface_tests.cpp` object only.
- A hand-linked focused executable passed 8/8 selected tests: empty vector/SHA
  boundaries, x and q reduction, platform Fp2 multiplication/squaring, x/q
  382-bit operations reduced against portable C++, and 64/128-bit
  quotient/remainder.
- Win64 portable assembly objects passed their exact 38-entry fallback-symbol
  allowlist; non-portable objects exposed zero fallback `$1` entries.
- The declaration inventory is exactly 77.
- Structural checks find no provider-renaming macros, `_MSC_VER` typo, stale
  commented declarations, or temporary ABI annotation leak.
- Native and portable MSVC `/W4 /WX` selector probes compile. Fixed native
  vector/SHA calls are direct tail calls.
- No full suite, sanitizer sweep, portable-AVX build, CUDA build, BTx/BTE test,
  or benchmark was run for this packet.

## Automated evidence

| Tool/profile | Result | Blind spots |
| --- | --- | --- |
| MSVC 19.50 compile/codegen | Strict native and portable selector probes pass; fixed native counts tail-call raw assembly | Does not validate Apple object format or AArch64 assembler linkage |
| Win64 symbol allowlist | Portable 38/38 exact entries; non-portable zero fallback entries | Public symbols outside the `$1` fallback relationship are covered by focused linking |
| Focused GoogleTest object/executable | 8/8 selected assembly-interface tests passed | q public dispatch can select x internally on ADX-capable hardware |
| CodeQL baseline | No declaration-header result beyond linked findings | Does not reason about handwritten assembly loop preconditions or timing |

## Residual risk and decision

The provider declaration and selection review is complete at
`2333650dcd3017f34e169c79bc8dc86e27d7305e`. The declaration layer now states
the operand ranges, modulus, alias rules, shift-count domain, and division
estimate domain. The raw positive-length prefetch symbol is
`vec_prefetch_asm`; the selected wrapper handles zero before assembly.

The exact MSVC object exports the new symbol and omits the obsolete raw name.
Independent x86-64 arithmetic tests pass 6/6. Apple ARM64 execution remains
required for the changed shared ARM prefetch owner. Non-Apple ARM64 and CHERI
remain experimental, and CUDA is outside this review. The manifest row remains
`in_review` until the generated assembly packet receives that evidence.

## Batch C Montgomery boundary follow-up

The declaration boundary was reread with all six Montgomery generators and
their selected callers at `a29895f`. It now states the sparse 256-bit operand
asymmetry, REDC and from-Montgomery bounds, exact alias rules, positive public
repeated-square count, Fp2 component order, 382x output shape, and timing
scope. These comments add no executable code.

The independent byte oracle passes for the selected, forced-generic, and
direct ADX providers on MSVC x86-64. `AUD-370` is validated. The header remains
in review with `AUD-369` and `AUD-371` pending supported Apple ARM64 evidence.

## Batch D inverse boundary follow-up

The inverse and residue-test declarations were reread with all seven owners at
`83c11ba`. The inverse functions return a bounded double-width almost-inverse,
not a conventional inverse. For `R=2^w`, a nonzero invertible input satisfies
`ret * inp = R^2 (mod mod)`; zero returns zero. The residue test requires a
reduced input and odd prime modulus and includes zero among quadratic residues.

These declaration changes add no instructions. Independent byte-arithmetic
and Euler-criterion oracles pass for selected, portable, direct ADX, and
forced-generic x86 providers. `AUD-372` is validated. Apple execution of the
new oracle remains the production evidence gap recorded by `AUD-373`.

## Batch E SHA boundary follow-up

The SHA declarations were reread with both retained generator owners at
`2f59883`. The raw compression contract now requires positive consecutive
blocks, a count at most `SIZE_MAX/64`, and disjoint state and message storage.
Copy and emit helpers permit disjoint or exact-alias storage, but not partial
overlap. Lengths and CPU capabilities are public; message bytes do not select
control flow or addresses.

The internal portable-C provider is declared as the x86 no-ISA fallback.
Forced C, SSSE3, and available SHA-NI known-answer tests pass. These comments
and the fallback declaration add no hot-path conversion or allocation.
`AUD-377` is validated.
