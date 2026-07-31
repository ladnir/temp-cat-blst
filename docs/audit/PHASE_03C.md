# Phase 3C Accepted Assembly Chain Evidence

## Record

- Date: 2026-07-30
- Baseline commit: `498baac`
- Implementation commit: `5d31975`
- State: validated and committed
- Scope: upstream `525a87a`, `4a74ffd`, `e789f8a`, `e3d708b`,
  `33d7f80`, and `bf008ec`; local `AUD-017`
- Benchmarks: disabled and not run
- CUDA: disabled

## Upstream translation

The accepted commits were applied in dependency order to the repaired
generator trust chain:

1. Remove the unused specialized Montgomery `sqr_mont_382x` implementation.
2. Preserve the upstream ADX carry-chain schedules exactly.
3. Update the x86 translator to retain braced directive tokens and correct
   the CFI spelling in its example.
4. Clarify the sparse 256-bit Montgomery operand contract.
5. Put the reduced, sparse `R^2` constant first in native C++ conversion and
   reciprocal calls.
6. Regenerate the complete output tree once and accept only explained files.

The ordinary non-Montgomery `sqr_382x` implementation remains in production
and is still used by the final exponentiation path.

## Source and generated-output provenance

After normalizing line endings, the five affected upstream-derived generator
or translator sources match upstream `41b511a` exactly:

- `blst/asm/mul_mont_384-armv8.pl`
- `blst/asm/mulq_mont_384-x86_64.pl`
- `blst/asm/mulx_mont_256-x86_64.pl`
- `blst/asm/mulx_mont_384-x86_64.pl`
- `blst/asm/x86_64-xlate.pl`

Regeneration changed exactly 20 checked-in outputs. Seventeen contain the
accepted routine removal or carry-chain changes. Three additional Win64
outputs contain only the Phase 3B guard-whitespace normalization required to
match the upstream generator spelling. After normalizing line endings, all
20 match upstream `14cd09a` exactly.

`build/verify-assembly.sh` generated two independent trees and verified all
100 outputs against the checked-in allowlist after the port.

## Fork-specific defect found during integration

Instantiating the native `blst_256_t::to()` and `reciprocal()` paths exposed a
latent fork regression recorded as `AUD-017`. The fork's assembly header had
moved generic 256- and 384-bit declarations into the ARM-only block and had
lost the upstream-style compile-time mapping from generic names to ADX names.
It also misspelled the MSVC predefined macro as `__MSC_VER`.

Consequences were:

- supported native C++ template members failed to compile when instantiated;
- non-portable ADX builds could not select direct ADX entry points through
  the generic native C++ surface; and
- portable builds lacked the declaration layer needed to enter q-function
  runtime dispatch.

The header now declares the generic surface for every assembly architecture.
Non-portable ADX builds map it directly to x entry points; portable and
non-ADX builds retain generic q entry points. The MSVC declaration guards now
use `_MSC_VER`.

## Validation

Toolchain:

- CMake with Ninja;
- MSVC C/C++ 19.50.35728.0;
- MASM 14.50.35728.0;
- `dumpbin` from MSVC 14.50.35717.

| Profile | Configuration | Result |
| --- | --- | --- |
| Scalar direct ADX | Debug, `BLST_PORTABLE=OFF`, `BLST_ADX=ON`, assembly on, AVX-512 off | 268 passed, nine AVX-only tests skipped, zero failed (277 registered) |
| Scalar portable dispatch | Debug, `BLST_PORTABLE=ON`, `BLST_ADX=OFF`, assembly on, AVX-512 off | 310 passed, nine AVX-only tests skipped, zero failed (319 registered) |
| AVX-512 portable | Debug, `BLST_PORTABLE=ON`, scalar default off, portable AVX-512 default on | 142/142 tests passed |

New native C++ regressions exercise unreduced boundary values through
`blst_256_t::to()` and compare reciprocal results with the scalar field
implementation.

The Win64 object allowlist contains exactly 38 fallback entries in portable
mode and none in non-portable mode across all 13 MASM objects. The two removed
entries were the consumer/provider pair for `sqr_mont_382x`.

Source search and `dumpbin /symbols` confirmed that `sqr_mont_382x` is absent
from all three validated libraries, while `sqr_382x` remains present.
`git diff --check` passed.

## Residual validation

- Win64 x86-64 scalar and AVX-512 code was assembled and executed.
- ARM64 and CHERI outputs match the accepted upstream generated artifacts but
  were not assembled or executed on this host.
- The tests are Debug builds; the generator schedules and assembly symbol
  surfaces are independent of C/C++ optimization level.
- Allocation behavior and performance were not changed or benchmarked in this
  phase.
