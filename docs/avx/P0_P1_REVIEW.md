# AVX Proof P0/P1 Closure Review

## Result

The four foundational correspondence rows are accepted after a structured
source review and focused executable checks. Acceptance means the exact C++
sources and generated constants refine the proved models at their recorded
hashes. It does not mean Lean proves Clang, the x86 ISA, object-file loading,
or the processor.

The review was performed on 2026-08-09. The portable semantics repair is
commit `63c57c63116051deabdbb2b89320853cd5e20d34`. Every mapped file hash is
enforced by `scripts/audit/validate_avx_correspondence.py`.

## P0 lane and IFMA mapping

The hand review followed values from the C++ wrappers through the portable
oracle and into the Lean definitions.

| C++ operation | Lean owner | Reviewed relation |
| --- | --- | --- |
| `vv512i` in native builds | `CatBLST.AVX.Vec8` | Exactly eight ordered 64-bit lanes; native `_vv512_*` wrappers call `_mm512_*` directly and contain no `cc` conversion |
| `_vv512_add_epi64`, `_vv512_sub_epi64` | `Word64.add`, `Word64.sub` | Independent lane-wise arithmetic modulo `2^64` |
| logical shifts and masks | `Word64.shiftLeft`, `Word64.shiftRight`, `Word64.bitAnd` | Lane-local logical operations; fixed production counts are in range |
| `_vv512_permutexvar_epi64` | `Vec8.permute` | Destination lane `i` reads the source lane selected by the low three bits of index lane `i` |
| lane-zero extraction and broadcast | `Vec8.lowLane`, `Vec8.broadcastLow` | The low 64-bit lane is extracted and copied to all eight lanes |
| `_vv512_madd52lo_epu64` | `madd52loLane` | Low 52 bits of both operands are multiplied; low product digit is added modulo `2^64` |
| `_vv512_madd52hi_epu64` | `madd52hiLane` | Low 52 bits of both operands are multiplied; high product digit is added modulo `2^64` |
| portable 26-bit split | `PortableIFMA.sourceLow`, `PortableIFMA.sourceHigh` | C++ `p0`, `p1`, `p2`, carry, mask, shift, and OR schedule matches the Lean transcription; every non-accumulator temporary fits `uint64_t` |

The portable C++ schedule uses unsigned arithmetic for every product,
temporary, and accumulator. `p0` and `p2` are below `2^52`; `p1` is below
`2^53`; the assembled IFMA digit is below `2^52`; accumulator addition wraps
as required. The high bits of each multiplicand are deliberately ignored.

The review found `AUD-224`: immediate shifts by 64 or more invoked C++ shift
undefined behavior instead of the AVX saturating semantics, and the unaligned
load/store helpers accessed storage through `uint64_t*`. Commit `63c57c6`
defines all shift-count cases, implements arithmetic right shift with unsigned
operations, and uses fixed-size `memcpy` for byte-safe loads and stores. These
changes affect only portable emulation; native `_vv` production wrappers are
unchanged.

### P0 executable evidence

- Local MSVC 19.50 compiled the changed emulation and compatibility-test
  translation units. The old focused MSVC target still has an unrelated
  standalone scalar-assembly link omission, so it is not counted as execution
  evidence.
- Peach AMD Ryzen 9 7950X, Clang 21.1.8, Release, `-march=native`, native and
  portable AVX-512 enabled: `29/29` focused compatibility and MatrixNoK tests
  passed. This includes deterministic boundary and randomized IFMA cases,
  masked high inputs, accumulator wrap, all immediate-shift boundaries through
  255, truly unaligned byte addresses, masks, comparisons, permutations,
  extraction, and broadcast.
- Clang 21 optimized lane-zero extraction to one `movq` load and return.
  Quotient broadcast compiled to `vpsrlq`, `vpaddq`, `vpbroadcastq`, and
  `vmovdqu64`, followed only by `vzeroupper` and return. There are no wrapper
  calls, conversion objects, or memory round trips.

## P1 constants and regeneration mapping

The generated header parser was reviewed as a deliberately narrow parser, not
as a C++ evaluator:

- every accepted scalar and array name is explicitly enumerated;
- the generated header contains decimal, nonnegative integer literals only;
- each initializer is delimited before numeric extraction;
- dimensions and 64-byte alignment declarations are checked independently;
- missing, renamed, truncated, reordered, or extra-shape data fails closed;
  and
- the generated Lean snapshot includes every parsed literal and the SHA-256 of
  the complete C++ header.

`validate_matrixnok_constants.py` independently reconstructs both RNS bases,
their pairwise-coprimality relations, CRT factors, Montgomery factors,
rotations, four conversion matrices, correction vectors, quotient estimates,
fixed encodings, and Fp12 offsets. Its corruption self-test independently
changes a product, rotation, matrix entry, permutation entry, quotient
estimate, and Montgomery factor; all are rejected.

The matching Lean packet proves the fixed dimensions and ranges, exact basis
products, inverse witnesses, permutation relation, all four regeneration
families, fixed-point error bounds, and the 128-bit quotient-accumulator bound.
Tower constants are connected separately by the already accepted P3 Frobenius
and cyclotomic rows; they are not silently claimed by the base-field P1
theorems.

### P1 executable evidence

- `validate_matrixnok_constants.py --self-test` passed.
- `generate_matrixnok_lean.py --check` reproduced the checked-in Lean snapshot
  byte for byte and confirmed source SHA-256
  `1643256db5038b0a643ded8053dc833ab2613d20dbc6de8122296f0a4699361a`.
- The focused Lean build of `Lane`, `IFMA`, `PortableIFMA`, `Constants`, and
  `Regeneration` completed successfully with the pinned toolchain.
- Native execution is not required to validate finite constant regeneration;
  native consumers remain covered by the accepted P2-P4 differential rows.

## Accepted boundary

| Row | Decision | Remaining external boundary |
| --- | --- | --- |
| `AVX-P0-001` | Accepted | Clang/x86 intrinsic lowering and processor execution |
| `AVX-P0-002` | Accepted | C++ implementation correctness beyond the reviewed source hash |
| `AVX-P1-001` | Accepted | Compiler placement/loading of immutable integer objects |
| `AVX-P1-002` | Accepted | Correctness of future parser or generator changes beyond their reviewed hashes |

Any change to a mapped C++ file or constant snapshot invalidates its recorded
hash and reopens the affected row. Constant-time machine behavior remains a
separate code-generation and leakage-testing claim; it is not established by
this closure.

The `AVX-P0-001` row was reaccepted at `9e237456`. The selected carrier remains
the exact eight-lane construction modeled by Lean. Exact lane equality replaces
two signed comparison masks, while mapped extraction, broadcast, and IFMA
arithmetic are unchanged. Focused portable execution and optimized native
inspection confirm the retained behavior. `AVX-P0-002` is rebound to the same
commit after removal of zero-consumer emulation operations; its reviewed IFMA
decomposition is unchanged.
