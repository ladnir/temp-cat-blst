# Phase 4 Allocation-Free Scalar MSM Precompute Evidence

## Record

- Date: 2026-07-30
- Baseline commit: `762e9a5`
- Portable assembly correction: `bab906e`
- MSM implementation commit: `c585f98`
- State: validated and committed
- Scope: `AUD-005`, upstream `84926f2`, and integration finding `AUD-018`
- Benchmarks: disabled and not run
- CUDA: disabled

## Scalar MSM root cause

The fork had replaced upstream's potentially excessive stack allocation in
the fixed-table MSM precompute tail with `malloc`. That avoided unbounded
stack growth but introduced a hidden, potentially multi-megabyte allocation
behind a void API. Arithmetic overflow or allocation failure silently zeroed
the output table, which could make failure look like a valid all-identity
precompute result.

The existing fast path already uses unused output-table capacity as temporary
projective storage. Only the final rows required a bounded replacement.

## Allocation-free design

The upstream `84926f2` horizontal-row strategy was translated into the fork's
C++ field and point types while retaining its checked arithmetic:

1. Validate the window, pointers, shifted table count, and output byte count
   before writing.
2. Preserve the existing in-table batch conversion while sufficient output
   capacity remains.
3. If all remaining projective rows plus the inversion accumulator fit under
   `SCRATCH_LIMIT`, allocate exactly that bounded amount on the stack and use
   the original whole-tail conversion.
4. Otherwise allocate at most `SCRATCH_LIMIT` bytes for one even-sized
   projective row chunk.
5. Convert the first chunk to affine in the output table, extend subsequent
   chunks from the last projective multiple and the first affine chunk, and
   batch-convert each chunk directly into its final position.

The default stack bound is 144 KiB, or 45 KiB for WebAssembly. A compile-time
check requires it to hold at least two projective points. The whole-tail size
also retains the fork's additional one-field accumulator padding rather than
copying the upstream allocation expression without that safety margin.

The output layout and the supported convention of storing packed input points
at the end of the output table are unchanged.

## Portable assembly integration finding

The Linux ASan validation exposed `AUD-018`. In portable mode,
`blst/backend/backend_asm.h` correctly selected generic q entry points, but
`blst/generated/assembly.S` included only the x/ADX providers. Linux therefore failed
to link native C++ field template calls such as `mul_mont_sparse_256`.

The wrapper now follows the upstream portable inclusion policy for ELF, COFF,
and Mach-O:

- direct ADX mode includes the x providers;
- non-ADX mode includes the q implementations; and
- portable mode includes both q dispatchers and x providers.

The validated ELF object exports generic q names and x names together with
the guarded `$1` fallback providers. This correction is isolated in
`bab906e`.

## Validation

Focused precompute coverage uses three inputs including affine infinity and
compares every table entry with an independently batch-normalized projective
row. It then evaluates the table and compares the MSM result with direct
scalar multiplication.

Covered cases:

- G1/Fp and G2/Fp2;
- packed, pointer-array, and output-table-aliased packed inputs;
- below-threshold whole-tail conversion;
- horizontal chunking above the threshold;
- multiple horizontal chunks at the maximum supported `wbits=14`;
- null output and overflowed size queries; and
- direct-ADX, portable-dispatch, AVX-512-default, Release, and Linux ASan
  configurations.

MSVC Debug's CRT allocation hook surrounded only the production precompute
call and observed zero allocations in every case. `dumpbin` and `nm` also
found no `malloc`, `free`, or C++ allocation references in the Debug, Release,
or GCC/ASan `multi_scalar` objects.

| Profile | Configuration | Result |
| --- | --- | --- |
| Windows scalar direct ADX | Debug, `BLST_PORTABLE=OFF`, `BLST_ADX=ON` | 270 passed, nine AVX-only tests skipped, zero failed (279 registered) |
| Windows scalar portable dispatch | Debug, `BLST_PORTABLE=ON`, `BLST_ADX=OFF` | 312 passed, nine AVX-only tests skipped, zero failed (321 registered) |
| Windows AVX-512 portable | Debug, portable AVX-512 default | 144/144 tests passed |
| Linux scalar portable dispatch | GCC 13 Debug, ASan and leak detection, q+x assembly | 313 passed, nine AVX-only tests skipped, zero failed (322 registered) |
| Windows optimized scalar | Release, direct ADX | Three focused threshold/invalid-input tests passed |

`scripts/build/verify-assembly.sh` verified all 100 deterministic generated outputs.
The audit manifest, upstream manifest, and `git diff --check` passed.

## Residual validation

- ELF portable assembly was built, symbol-inspected, and executed. The same
  wrapper correction for COFF non-MSVC and Mach-O mirrors upstream but was not
  built on this host.
- ARM64 and WebAssembly were not executed. Their scalar algorithm uses the
  same bounded logic; WebAssembly selects the smaller 45 KiB limit.
- This phase reviewed the fixed-table precompute section, its C API, and its
  tests. It does not mark the entire `multi_scalar.cpp` translation unit
  complete in the file-by-file audit.
- No performance benchmark was run because benchmarks are outside the audit
  scope.
