# File Review: `blst/backend/backend_portable.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/backend/backend_portable.h` |
| Roadmap phase / queue | Phase 13, portable-backend trust chain |
| Reviewed commit | `2333650dcd3017f34e169c79bc8dc86e27d7305e` |
| Reviewer | Codex audit team |
| Date | 2026-08-13 |
| Disposition | `in_scope_production` |
| Subsystem | 64-bit portable and ADX limb primitives, arithmetic declarations, constant-time helpers, SHA-256 fallback |
| Profiles | scalar C++, retained C/Rust consumers; experimental CUDA declarations are not release-gated |
| Generation/provenance | Handwritten fork file with upstream BLST arithmetic lineage |
| Related files | `backend_portable.cpp`, `backend.h`, `backend_asm.h`, scalar fields, `common/defines.h`, `test/Interface_tests.cpp` |

## Purpose and trust boundary

This header defines the lowest-level two-word arithmetic representation,
portable implementations of single-limb operations, ADX intrinsic mappings,
constant-time bit helpers, fixed-width arithmetic declarations, prefetching,
and portable SHA-256 primitives. Limb values can contain secret key and field
material. Supported host paths must be allocation-free, branch-independent of
secret data where the algorithm is intended to be constant-time, and directly
inlineable into arithmetic kernels.

## Contracts and invariants

- CatBLST uses 64-bit `limb_t` exclusively.
- `wide_limb_t` is exactly two contiguous 64-bit words in low/high order.
  Its private word type is `unsigned long long`, matching Intel intrinsic
  output-pointer types exactly on both LP64 and LLP64.
- The portable namespace is the pure C++ oracle. ADX functions map directly to
  `_mulx_u64`, `_addcarry_u64`, and `_subborrow_u64`; selection is compile-time.
- Carry and borrow inputs are zero or one. Subtraction high words encode borrow
  as zero or all ones where documented.
- `cc_num_bits(0)` is zero on every compiler. MSVC scans `l | 1` so the output
  parameter is initialized without branching on the input.
- Fixed-width template declarations inherit their array-size and alias
  preconditions from their definitions and callers.
- No supported host operation allocates, throws, performs runtime dispatch, or
  owns caller storage.

## Batch B contract delta

`AUD-368` makes the existing narrow contracts explicit. Modular add,
subtract, shift, small multiplication, and conditional negation require
reduced operands. Right shift also requires an odd modulus. Checked scalar
addition and subtraction consume canonical little-endian byte strings.
Division correction consumes an exact or one-high estimate. These changes are
documentation only and add no generated code.

## Review evidence by checklist area

### Identity, lineage, and dependency map

The file was reviewed with `backend_portable.cpp`, the scalar selector,
assembly declarations, field callers, and retained tests. The Clang 8-11
funnel-shift workaround remains narrowly versioned to upstream miscompile
`a061960b`; modern Clang and AppleClang keep the forced-inline path.

### Correctness and undefined behavior

The review covered multiply, add/carry, subtract/borrow, shifts, modular
helpers, inversion support, conditional masks, wide layout, and SHA fallback
interfaces. `AUD-052` removed incompatible output-pointer casts rather than
relying on equal width. `AUD-054` removed the uninitialized MSVC bit-scan
output. The last unreachable 32-bit implementation and selector were removed
under `AUD-047`. Size, offset, width, and alignment assertions protect the
wide representation.

### Memory, lifetime, API, and concurrency

The header owns no dynamic storage and retains no pointer. All state is caller
provided or automatic. It exposes internal C++ declarations, not a supported
public ABI; nevertheless its fixed representations feed the C and Rust ABI.
There is no mutable global, lazy, cached, or thread-local state.

### Side channels and secret lifecycle

Carry/borrow and mask helpers remain fixed instruction sequences. The MSVC
zero-bit scan compiles to `or; bsr; test; cmovne`, with no input-dependent
branch, address, call, or spill. The ADX wrappers compile directly to arithmetic
instructions and caller-output stores. Unused diagnostic printing and the
header's transitive stdio dependency were removed. This file neither creates
entropy nor persists secrets.

### Portability and diagnostics

Exact post-cleanup objects compile under MSVC `/W4 /WX` and GCC/Clang
`-Wall -Wextra -Wpedantic -Werror`. The translation-unit-wide Clang diagnostic
suppression was removed. CUDA remains experimental and was not built or used
as release evidence in this packet.

### Performance-sensitive structure

Operations remain forced-inline, fixed-width, stack-only, and selected at
compile time. No generic runtime abstraction, conversion temporary, memory
round-trip, or dispatch branch was introduced. MSVC Release wrappers contain
direct `mulx`, `adc`, and `sbb`; a benchmark was unnecessary because the
instruction mapping was inspected directly.

## Findings

- `AUD-036` - validated: modern Clang inlining/warning behavior.
- `AUD-047` - validated: unsupported 32-bit branches removed.
- `AUD-052` - validated: exact ADX output typing.
- `AUD-054` - validated: defined, branchless MSVC zero bit scan.

## Validation

- 72 deterministic boundary combinations plus 10,000 randomized combinations
  compare all ADX wide primitives with the pure portable implementation under
  MSVC 19.50, GCC, and Clang.
- MSVC Release assembly for multiply/add/subtract contains direct arithmetic
  and output stores only; the `cc_num_bits` probe has no secret-dependent
  branch, call, address, or spill.
- Focused portable inverse output matches assembly for zero, one, and
  modulus-minus-one.
- Exact backend objects pass strict MSVC, GCC, and Clang compilation.
- `git diff --check` and the audit manifest validator pass at closure.

## Residual risk and decision

Complete for the supported host contracts at the reviewed commit. Constant-
time and code-generation evidence is compiler/version/option-specific and must
be refreshed for supported toolchain changes. Experimental CUDA behavior is
not promoted by this review. Broader field/protocol correctness, sanitizer,
fuzzer, and full release matrices remain their subsystem gates.

## Batch D inverse follow-up

At `83c11ba`, the declarations now distinguish the double-width
almost-inverse from a conventional inverse. They state the reduced input,
odd-modulus, left-aligned modulus, disjoint-array, zero, bound, congruence, and
fixed-schedule contracts. The residue test states its odd-prime domain and
includes zero among residues. These documentation changes add no instructions.
