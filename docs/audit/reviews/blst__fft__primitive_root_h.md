# File Review: `blst/fft/primitive_root.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/fft/primitive_root.h`; compatibility shim `blst/fft/primitveRoot.h` |
| Roadmap phase / queue | Phase 17 FFT support-header packet |
| Reviewed commit | `4a6e95126f6977e99b2330a6753fa66fab1034b2` |
| Reviewer | Codex audit team |
| Date | 2026-08-11 |
| Disposition | `in_scope_production` |
| Subsystem | Deterministic power-of-two FFT root construction and validation |
| Profiles | scalar and selected AVX native C++; generic explicit-generator helper |
| Generation/provenance | Locally maintained CatBLST header with BLS12-381 constants independently regenerated from field parameters; Apache-2.0 and Category Labs |
| Upstream references | `not_applicable: upstream BLST has no corresponding FFT root-table helper` |
| Related files | `blst/fft/domain.h`, `blst/fft/ctFft.h`, `blst/fft_c_api.cpp`, `scripts/audit/validate_fft_roots.py`, `test/Fft_tests.cpp`, `test/Fft_api_tests.cpp` |

## Purpose and trust boundary

The authoritative header supplies exact power-of-two roots for the BLS12-381
scalar field, validates exact root order, optionally derives a root from an
explicit caller-supplied field generator, and prepares ordered power tables.
Transform sizes and roots are public domain parameters. The helper handles no
secret entropy and performs no allocation.

The historical misspelling is retained only as a one-include compatibility
shim for frozen BTx/BTE and excluded experimental CUDA source. Supported
CatBLST sources include the correctly spelled header.

## Contracts and invariants

- Built-in Fr roots exist only for nonzero powers of two through `2^32`.
- `isPrimRootOfUnity` returns false for invalid order shapes, zero, and roots
  whose order is a proper divisor of the requested power of two.
- Other fields must supply a generator explicitly. The order must divide the
  field multiplicative-group order, and the derived root is checked exactly.
- `precomputeRootsOfUnity` validates its nonempty power-of-two span and exact
  root order before modifying output. Invalid input throws
  `std::invalid_argument` with no partial publication.
- `psi` and output storage do not overlap. Valid preparation writes
  `1, psi, ..., psi^(n-1)` in order.

## Dependency and call-path map

`domain.h` and `middle_product.h` consume `primRootOfUnity<Fr>` for supported
public sizes. Ordinary FFT tests and AVX correspondence tests inspect the same
table. The generic preparation helper is retained for boundary tests and
caller-supplied finite-field generators. The C FFT boundary catches native
exceptions and reports failure. Frozen BTx/BTE compiled through the unchanged
misspelled include and was not audited or modified.

## Checklist evidence

### Identity, lineage, and API

The full history and every tracked caller were inventoried. The new header has
standalone dependencies, Apache-2.0/Category Labs provenance, deterministic
APIs, `[[nodiscard]]` query/build results, and explicit invalid-input behavior.
The old path exports no independent declarations or state.

### Correctness

For a power-of-two order, exactness is established by `u^n = 1` and
`u^(n/2) != 1`; the Fr specialization uses one repeated-squaring chain. The
independent Python validator derives `7^((r-1)/2^32)`, repeatedly squares in
the canonical field, Montgomery-encodes each result, and matches all four
little-endian limbs of all 33 entries. Its corruption self-test fails closed.

### Memory, lifetime, and undefined behavior

All state is immutable static storage or caller-owned span storage. Size and
order validation precede indexing and writes. No zero shift, factor-list
underflow, random-search nontermination, dynamic initialization, or partial
output remains. Native overlap is a documented precondition.

### Side channels, randomness, and concurrency

Sizes, exponents, indices, and roots are public FFT parameters. There is no
random discovery or secret lifecycle. The root table is constant-initialized
immutable storage and can be shared without synchronization.

### Portability and build behavior

The implementation uses C++20 bit operations already required by the project,
fixed-width limbs, and no host-endian exponent byte reinterpretation. The
correct header compiles standalone under MSVC 19.50 at `/W4`. The compatibility
shim successfully rebuilt its frozen consumers.

### Performance-sensitive structure

The common valid Fr lookup has no field arithmetic or dynamic initialization.
Optimized MSVC Release code performs public shape/range checks, `tzcnt`, and one
indexed 32-byte load/copy. Exact domain validation replaces the former Release
root exponentiation with a bounded squaring chain and remains outside FFT
butterfly loops. No benchmark was warranted because the lookup strictly removes
first-call work and its object code is directly verified.

## Findings

`AUD-043`, `AUD-049`, and `AUD-306` through `AUD-308`.

## Validation

- Independent root regeneration plus deliberate-corruption self-test: passed
  all 33 entries.
- MSVC 19.50 standalone C++20 `/W4` compile: passed without diagnostics.
- Scalar Debug and Release focused root/API tests: 5/5 passed in each profile.
- Optimized Release symbol/disassembly inspection: immutable table present;
  legacy `$TSS`/root lambda/guard absent; valid lookup has no calls.
- Frozen BTx/BTE compiled unchanged as a dependency. Its existing unrelated
  ignored-`nodiscard` warning remains.
- No benchmark, CUDA run, full suite, sanitizer, native AVX-512, remote, or
  Apple run was performed.

## Residual risk and decision

The table validator is source-level independent evidence; it does not replace
compiler, hardware, sanitizer, fuzzing, native AVX-512, or Apple milestone
gates. The generic helper assumes its field type correctly implements the
documented finite-field operations. This owner and its compatibility shim are
complete for line-by-line hand review.
