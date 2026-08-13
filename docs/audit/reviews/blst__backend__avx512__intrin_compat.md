# File Review: AVX-512 intrinsic compatibility boundary

## Record

| Field | Value |
| --- | --- |
| Path | `blst/backend/avx512/intrin_compat.h`; `blst/backend/avx512/intrin_compat.cpp` |
| Roadmap phase / queue | Phase 16 selected-vector boundary follow-up |
| Reviewed commit | `9e237456a946abe827aa55a8a122e4a040114ffc` |
| Reviewer | Codex audit team |
| Date | 2026-08-13 |
| Disposition | `in_scope_production` |
| Subsystem | Selected native AVX-512 operations and developer-only C/C++ emulation |
| Profiles | native AVX-512/IFMA; portable AVX Debug |
| Generation/provenance | Hand-maintained Category Labs code introduced during the CatBLST AVX integration; operation semantics follow Intel AVX-512 intrinsics |
| Upstream references | `not_applicable: upstream BLST has no selected-vector compatibility owner`; local history starts at `52cfdf66` and the native/portable split at `0b7ec477` |
| Related files | `blst/avx/cpu/vector/avx512ifma.hpp`, `test/avx512_compat_test.cpp`, MatrixNoK field and curve owners |

## Purpose and trust boundary

The header defines the selected `vv` vector domain. A native build maps each
operation directly to an Intel intrinsic. A portable AVX Debug build maps it
to an independent pure C/C++ operation with `cc` storage. A combined build
keeps both domains but selects native vectors for production calls.

The implementation file owns only the C/C++ emulation. This vector emulation
is a developer oracle and fallback for Debug testing. It is not a deployment
backend or evidence about native constant-time code generation. This does not
change scalar terminology: scalar `cc` remains the supported production
C/C++ implementation.

Lane values can depend on secrets in native production calls. The selected
wrappers must therefore add no dispatch, conversion, allocation, exception,
or memory round trip.

## Contracts and invariants

- `cc512i`, `cc256i`, and `cc128i` exist only with
  `BLST_AVX512_PORTABLE`. They are standard-layout, trivially copyable arrays
  of eight, four, and two `uint64_t` lanes.
- `vv512i`, `vv256i`, and `vv128i` alias native intrinsic types whenever
  `BLST_AVX512` is enabled. Native selection wins in combined builds.
- Native execution requires an OS context and CPU supporting AVX-512F,
  AVX-512DQ, and AVX-512IFMA.
- Unaligned loads and stores require 64 readable or writable bytes. Masks use
  their low eight bits. Variable permutation indices use their low three bits.
- Immediate extraction indices are compile-time checked. Immediate shifts
  outside `[0,63]` produce zero, matching the retained Intel operations.
- IFMA operations truncate each multiplicand to 52 bits and add the selected
  half-product modulo `2^64`.
- Equality compares unsigned lane bit patterns and returns an eight-bit mask.
- All operations are allocation-free, nonthrowing, and use no shared mutable
  state. Exact input/output aliasing is harmless because operations take
  values or complete the input read before a store.
- Native-to-`cc` bridges exist only in combined-build differential tests.
  Production selected-vector code must not call them.

The files expose no public C or Rust declaration and promise no stable ABI.

## Dependency and call-path map

MatrixNoK vector, RNS, field, curve, GT, pairing, FFT, and middle-product
owners consume the `vv` operations through `avx512ifma.hpp`. The portable
implementation is compiled only when `BLST_AVX512_PORTABLE` is enabled. The
combined compatibility target calls both domains directly and compares every
retained emulated operation with its native intrinsic.

Repository-wide call inventory at the reviewed commit found 16 compatibility
operations with no production consumer. Those operations and their tests were
removed. The retained operation set is exactly the set used by production or
needed at an intentional test bridge.

## Upstream/specification comparison

Upstream BLST has no corresponding selected-vector layer. The native behavior
is specified by the Intel intrinsic interfaces. The portable implementation
is an independent source-level oracle for the retained subset. It deliberately
does not emulate unrelated AVX-512 operations.

The deleted `avx_layout.h` belonged to CatBLST's former packed base-field
implementation. Commit `75df5873` removed that implementation but left the
562-line, zero-caller helper layer behind. No live owner included the file at
the reviewed commit.

## Checklist evidence

### 1. Identity and purpose

Both files carry Apache-2.0 and Category Labs ownership. Their native and
developer-only roles are explicit.

### 2. Contracts and invariants

The representation, selection, CPU, byte-count, mask, index, shift, IFMA,
alias, allocation, exception, and state contracts are stated above.

### 3. Upstream and specification lineage

The owner is CatBLST-specific. Local introduction and separation commits are
recorded, and the retained operations follow Intel intrinsic semantics.

### 4. Correctness

Ten combined-build tests compare retained portable operations with native
intrinsics over deterministic boundary and random inputs. Dedicated IFMA
cases cover masked multiplicands, accumulator wrap, carries, and random lanes.

### 5. Memory and lifetime safety

The lane structures own their storage. Unaligned transfers use `memcpy` and
copy exactly 64 bytes. Compile-time layout assertions reject representation
drift.

### 6. Undefined and implementation-defined behavior

Unsigned lanes remove the former signed-conversion dependency from portable
comparison. Shift guards prevent C++ shifts by negative or excessive counts.
Native unsigned literals use `bit_cast` into the intrinsic's signed parameter
type, and the equal-width assumption is statically checked.

### 7. Side channels and cryptographic hygiene

Native `vv` wrappers are forced inline direct intrinsic calls. Peach GCC 15.2
emits one `vpcmpq` for selected equality and one `vpmadd52luq` for selected
low IFMA, byte-for-byte matching the direct-intrinsic probes. Neither probe
contains a call, spill, conversion, or memory access. Portable Debug is not a
constant-time claim.

### 8. Randomness and secret lifecycle

N/A: these fixed-width operations create no randomness and own no persistent
secret buffer.

### 9. API and ABI

N/A: the names and types are internal C++ implementation details.

### 10. Concurrency

Every operation uses values and caller-owned memory. There is no mutable
global state or lazy initialization.

### 11. Resource exhaustion and denial of service

Work and storage are fixed. No operation allocates or accepts an unbounded
count.

### 12. Portability and build behavior

Portable-only code compiles cleanly with Clang 18 and GCC 13 under
`-Wall -Wextra -Wpedantic`. Combined native-plus-portable execution passed on
Peach with GCC 15.2 and AVX-512/IFMA hardware. Portable AVX remains a hidden
developer Debug option.

### 13. Performance-sensitive structure

Native wrappers contain direct intrinsics only. Replacing two signed compare
operations with equality removes one `vpcmpq` from GCC 15 output. Removing the
dead layout and unused compatibility surface adds no production work.

### 14. Tests and evidence

The exact focused commands and exclusions appear below.

### Specialized additions

The reviewed header SHA-256 is
`867a8deefd15c95a13cb0df337ddfe66bfac26f506a78a1f386e9f98e00de434`.
The implementation SHA-256 is
`7fbe7d9c50ffeff8700b04d4e3254800e6fe09cd2100b850a87e65aafa1bbeb1`.

## Findings

`AUD-036`, `AUD-037`, `AUD-220`, `AUD-224`, `AUD-384`, `AUD-385`, and
`AUD-386`.

## Validation

- Structural search found no include of `avx_layout.h`, no production use of
  a native-to-`cc` bridge, and no retained spelling for any removed operation.
- WSL Clang 18 and GCC 13 compiled `intrin_compat.cpp` as C++23 portable-only
  code with `-Wall -Wextra -Wpedantic` and no diagnostic.
- Local MSVC 19.50 compiled the portable implementation and MatrixNoK
  primitive test objects, plus combined implementation and compatibility test
  objects, one target at a time.
- Peach GCC 15.2 Release built only `blst_avx512_compat_tests` with native and
  portable support enabled. All 10 discovered tests passed.
- Peach GCC 15.2 `-O3 -DNDEBUG -march=native` emitted identical selected and
  direct bodies for equality and low IFMA. Equality is one `vpcmpq`; IFMA is
  one `vpmadd52luq`.
- No benchmark, broad build, CUDA, BTx/BTE audit, portable AVX Release build,
  sanitizer, Apple CI, or full test suite ran.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing C/C++ baseline | `AUD-036`, `AUD-037` | Intrinsic lowering and constant-time machine behavior require direct inspection |
| Clang Static Analyzer / clang-tidy | Existing curated baseline | No additional retained-surface finding | AVX intrinsic semantics and generated code require differential tests and inspection |
| Clippy | N/A | N/A | No Rust source |
| ASan/UBSan | Not rerun | Deferred milestone evidence | Native and portable profiles require separate focused runs |
| Fuzzing | Not rerun | Deferred milestone evidence | The finite primitive surface has deterministic differential coverage |

## Residual risk and decision

The retained source and selected native lowering are complete for the reviewed
operation set. The portable oracle does not prove Intel hardware semantics,
compiler correctness, or constant-time execution on every supported compiler.
Those are explicit machine boundaries. Focused GCC 15.2 hardware execution and
code generation are sufficient to close this file review; broader release
profiles remain milestone evidence.

