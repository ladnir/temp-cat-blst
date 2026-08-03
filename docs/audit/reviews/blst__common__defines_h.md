# File Review: `blst/common/defines.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/common/defines.h` |
| Roadmap phase / queue | Phase 13, first configuration/representation packet |
| Reviewed commit | `80ae0f9a23da1c30415498241e204297e15157af` |
| Reviewer | Codex |
| Date | 2026-08-02 |
| Disposition | `in_scope_production` |
| Subsystem | common representation, compiler, assertion, and sizing contracts |
| Profiles | scalar, AVX-512, native C++ internals; C and Rust through compiled ABI |
| Generation/provenance | Locally composed C++ header with declared lineage from upstream BLST `src/vect.h` |
| Upstream references | Upstream head `54e6e55674722fc2797ebb4bbb71b26d881eb4b8`, `src/vect.h`; no one-to-one upstream file |
| Related files | `blst/config.h.in`, root `CMakeLists.txt`, portable/scalar/AVX backends, FFT helpers, `test/audit_launder_codegen.cpp`, `test/Fft_api_tests.cpp` |

## Purpose and trust boundary

This internal C++ header fixes CatBLST's fundamental word and vector
representations, compiler attributes, assertion/assumption behavior,
constant-time value barrier, host/device decoration, and two common sizing
helpers. Every arithmetic subsystem depends on its definitions. It must remain
zero-state and allocation-free, must not perform runtime backend dispatch, and
must not conceal compiler diagnostics across an including translation unit.

Inputs to `log2ceil` and `ceil_div` are public sizes. Limb/vector values may be
secret when used by downstream arithmetic, so representation and
`BLST_launder` must not introduce secret-dependent branches, addresses, calls,
or spills.

## Contracts and invariants

- `BLST_LIMB_SIZE` is exactly 64 and `limb_t` is exactly `uint64_t`.
- `uptr_t` is the standard `uintptr_t`, so supported 64-bit hosts have an
  integer type capable of representing object pointers.
- `NLIMBS(bits)` is used for widths divisible by 64; fixed vector aliases are
  contiguous limb arrays with no ownership.
- `BLST_ASSERT` terminates on failure in Debug and is absent in Release.
  `BLST_ASSUME` is therefore valid only where callers prove the condition in
  production; it is not an input-validation mechanism.
- `BLST_launder` is a compiler value barrier, not a memory fence or
  zeroization primitive.
- `log2ceil(0)` asserts in Debug and returns zero in Release; positive inputs
  return the least exponent whose power of two covers the input.
- `ceil_div(x,0)` asserts in Debug and returns fail-closed zero in Release;
  otherwise it computes the ceiling without add-before-divide overflow.
- `vec_zero` and `vec_is_equal` are declarations whose implementations and
  call sites are reviewed with the backend; this header owns no storage.

## Dependency and call-path map

The generated `config.h` supplies build choices. This header is included by
portable, assembly-selection, scalar-field, curve, signature, FFT, AVX-layout,
and API adapter code. Public C headers do not expose these C++ declarations,
but public C and Rust layouts depend on the compiled native representations
remaining compatible. FFT root preparation calls `log2ceil`; sizing and
scratch callers use `ceil_div`. Portable arithmetic and conditional operations
consume `BLST_launder`.

## Upstream/specification comparison

Upstream `src/vect.h` at the audited head defines analogous limb/vector aliases,
byte expansion, restrict compatibility, and GNU/Clang value laundering. This
header is not a wholesale copy: it is the local C++ representation boundary,
adds CMake feature consumption, host/device attributes, assertions,
fixed-width aliases, and sizing helpers. CatBLST intentionally rejects
upstream's dormant 32-bit/no-assembly limb representation. The Apache-2.0/SPDX
header and the local/upstream-lineage note now make that provenance explicit.

## Checklist evidence

### 1. Identity and purpose

The path, exact commit, Phase 13 position, local composition, upstream lineage,
supported profiles, and downstream consumers are recorded above. The header is
handwritten and has no generated output.

### 2. Contracts and invariants

The exact representation, assertion, assumption, barrier, and sizing contracts
are listed above. Zero and maximum sizing boundaries have focused tests. Debug
assertions supplement rather than replace the defined Release behavior of both
helpers.

### 3. Upstream and specification lineage

The file was compared with upstream head `src/vect.h`. Intentional divergences
are C++ namespacing/types, CMake-controlled feature attributes, 64-bit-only
limbs, local assertion policy, and common sizing helpers. The license header
identifies upstream lineage without representing the file as a wholesale copy.

### 4. Correctness

Fixed vector widths divide exactly by 64. `TO_LIMB_T` and `TO_BYTES` preserve
all eight bytes of a 64-bit word. `log2ceil` uses `std::bit_width(x - 1)` only
after the zero guard. `ceil_div` uses `x / y + (x % y != 0)`, avoiding the prior
addition overflow. No signed arithmetic is used in these size formulas.

### 5. Memory and lifetime safety

The header allocates and owns nothing. Vector aliases are caller-owned fixed
arrays. `uptr_t` is the standard integer pointer representation and is used
only for comparison/arithmetic preconditions. The assertions print static
metadata and terminate; they do not retain caller data. Downstream pointer and
array operations remain subject to their own file reviews.

### 6. Undefined and implementation-defined behavior

The non-64-bit preprocessor path is rejected. The former 64-bit shift after
`log2ceil(0)` is prevented by both the helper result and FFT boundary check.
Zero division is prevented before quotient/remainder evaluation. Compiler
attributes and inline assembly are selected explicitly per toolchain. Optimized
barrier behavior was inspected rather than inferred from source.

### 7. Side channels and cryptographic hygiene

`BLST_launder` probes for select, swap, conditional zero, secure zero, and a
portable arithmetic mask are branch/address/call/spill-clean for secret values
under MSVC 19.50, Clang 18.1.3, and GCC 13.3. The barrier itself was unchanged.
Sizing values and assertion paths are public. This file does not claim that all
downstream algorithms are constant-time.

### 8. Randomness and secret lifecycle

N/A: no entropy, nonce, secret creation, or persistent storage. The header
declares `vec_zero`; compiler-resistant clearing behavior is reviewed with its
implementation and callers.

### 9. API and ABI

The internal limb/vector layout is a compiled ABI dependency for C and Rust
opaque/fixed representations. Root CMake and this header independently enforce
64 bits. The file exports no public symbol and is unaffected by the future
build-time `blst_`/`catblst_` prefix except that both modes must share this exact
representation. No runtime prefix dispatch is present.

### 10. Concurrency

No global, static mutable, lazy, cached, or thread-local state exists. Assertion
output may interleave across failing threads, but failure terminates and is not
a production synchronization contract.

### 11. Resource exhaustion and denial of service

The helpers are constant time and allocate nothing. Invalid zero sizing returns
fail-closed zero in Release. Very large valid inputs do not overflow
`ceil_div`; `log2ceil` is bounded by the fixed-width bit operation. Callers must
still bound the work or allocation represented by returned sizes.

### 12. Portability and build behavior

The header compiled standalone and through focused backend/C API/FFT adapter
objects under strict MSVC, Clang, and GCC warning settings. Translation-unit-
wide Intel/GCC/MSVC diagnostic suppressions were removed, and exposed warnings
were repaired at their sources without a replacement pragma. CUDA remains an
experimental consumer; other unsupported configurations are outside the
production claim.

### 13. Performance-sensitive structure

Types remain fixed-width arrays and all small helpers are inline/constexpr.
The barrier probe found no added calls or spills. Replacing the GNU statement
expression with `cc_msb` retained forced inlining; strict optimized objects have
no helper call. No allocation, dynamic polymorphism, runtime dispatch, or hot
loop widening was introduced. No benchmark was warranted because optimized
structure was unchanged and directly inspected.

### 14. Tests and evidence

Focused helper/root tests cover zero, one, `SIZE_MAX`, Debug assertions, empty
roots, size one, and C failure translation. Strict compiler objects cover the
broad-warning-removal fallout. Exact compiler/hashes for value-barrier evidence
are in `PHASE_13.md`. Full suites and broad sanitizers are subsystem gates, not
required for this header-only packet.

### Specialized additions

N/A for assembly generation, FFT mathematics, and Rust unsafe implementation.
The AVX-512 relevance is type/configuration fan-out only; direct operation
selection was separately validated by `AUD-037`.

## Findings

`AUD-043`, `AUD-044`, `AUD-045`, `AUD-046`, and `AUD-047` - validated.

## Validation

- Standalone include compiled with MSVC `/W4 /WX` and Clang/GCC
  `-Wall -Wextra -Wpedantic -Werror`.
- Exact `backend_portable.cpp`, `c_api.cpp`, and `fft_c_api.cpp` objects compiled
  sequentially in each strict-warning profile.
- Focused Release and Debug helper/root tests passed 2/2 in each selected
  profile; the empty prepared-domain C API test passed 1/1.
- Explicit 64-bit CMake configuration succeeded; explicit 32-bit configuration
  failed before compilation.
- Optimized value-barrier assembly was inspected and hashed as recorded in
  `PHASE_13.md`.
- `scripts/audit/check_audit_manifest.ps1`: passed at closure.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Phase 12 scalar and portable-AVX databases | `AUD-046` source-hygiene context; no unresolved finding in this header | Constant-time and macro configurations require compiler evidence |
| Clang Static Analyzer / clang-tidy | Phase 12 retained profiles plus strict Clang 18 compile | No unresolved diagnostic | Cannot prove constant-time or downstream caller preconditions |
| Clippy | Retained Rust feature baselines | N/A for C++ header | Rust sees the compiled C ABI rather than this header |
| ASan/UBSan | Prior Phase 12 scalar baselines | No report mapped to these helpers | Sanitizers do not prove constant-time and Debug assertions intentionally terminate |
| Fuzzing | Deferred to subsystem expansion | N/A for finite helper boundaries in this packet | Future fuzzers should exercise downstream size composition and C adapters |

## Residual risk and decision

Complete for this header's contracts at the reviewed commit. The value-barrier
evidence is compiler/version/option-specific and must be refreshed when a
supported toolchain changes code generation. `BLST_ASSUME` remains a deliberate
production precondition primitive; every use must be justified during its
consumer's file review. Two 32-bit branches outside this file remain deferred
to `backend_portable.h`, and broader arithmetic/caller safety is not implied by
completion of the common definitions.
