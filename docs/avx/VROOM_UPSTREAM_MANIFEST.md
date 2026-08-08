# VROOM Upstream Manifest

## Baseline

- Repository: <https://github.com/SimonLangowski/VROOM>
- CatBLST import baseline: `31e98bbf7f32e1345299d4276dd6a5bc44dc0b69`
- Reviewed upstream head: `26bd9aee2168f167583cf096add355547a7b56d4`
- Delta: 31 commits
- Algorithm anchor: `d6ecd534c54a41702ea80ae75b52250b6d51b372`

Only the algorithm anchor changes executable CPU arithmetic. Later core-path
touches in `7995e57` and `8da9352` change comments only. The selected action is
a semantic BLS12-381 translation, not a merge.

## Decisions

| # | Commit | Summary | Core impact | Decision | State |
| ---: | --- | --- | --- | --- | --- |
| 1 | `50cb76d` | add BLST benchmarks | Benchmarks only | reject | closed |
| 2 | `88d58c0` | copy GPU code | GPU only | reject | closed |
| 3 | `d6ecd53` | BN254 plus camera-ready CPU results | MatrixNoK, generic RNS/curve changes, BN254, tests and benchmarks | translate BLS12-381 subset | in review |
| 4 | `49e2c8d` | add examples | Examples only | reject | closed |
| 5 | `7995e57` | artifact cleanup | `bounded_ring.hpp` comment cleanup only | already equivalent / document locally | closed |
| 6 | `07bc667` | reproducibility badge work | Artifact documentation | reject | closed |
| 7 | `31cc0af` | artifact test preparation | Artifact scripts | reject | closed |
| 8 | `f206c11` | update libraries | Artifact dependencies | reject | closed |
| 9 | `8da9352` | remove old benchmarks | Benchmark removal plus one comment | reject | closed |
| 10 | `b258f92` | simplify example | Examples only | reject | closed |
| 11 | `c177d97` | simplify example | Examples only | reject | closed |
| 12 | `0bf5966` | add examples | Examples only | reject | closed |
| 13 | `c35c0e5` | add examples | Examples only | reject | closed |
| 14 | `340c212` | add examples | Examples only | reject | closed |
| 15 | `20335ed` | fix artifact | Artifact only | reject | closed |
| 16 | `8a6c6d0` | artifact run | Results/artifact only | reject | closed |
| 17 | `2abf635` | warning and resource-usage cleanup | Artifact/tooling only | reject | closed |
| 18 | `3908993` | executable artifact files | Artifact only | reject | closed |
| 19 | `e0e7dcd` | remove dalek baseline | Baselines only | reject | closed |
| 20 | `4761fe0` | merge VROOM-test | Merge of artifact work | reject | closed |
| 21 | `58c8b23` | add missing artifact files | Artifact only | reject | closed |
| 22 | `6674d00` | specify dependency commits | Artifact metadata | reject | closed |
| 23 | `8619fe3` | GPU artifact pass | GPU/artifact only | reject | closed |
| 24 | `368c9e4` | merge VROOM-test | Merge of artifact work | reject | closed |
| 25 | `b8d936e` | metal-machine artifact run | Results/scripts | reject | closed |
| 26 | `fd6d149` | merge VROOM-test | Merge of artifact work | reject | closed |
| 27 | `f12e58d` | camera-ready README results | Documentation/results | use as performance motivation only | closed |
| 28 | `f2af184` | fix artifact toolchain | Tooling only | reject | closed |
| 29 | `45ed40a` | fix artifact toolchain | Tooling only | reject | closed |
| 30 | `b450276` | fix artifact toolchain | Tooling only | reject | closed |
| 31 | `26bd9ae` | artifact-evaluation feedback | Documentation/tooling only | reject | closed |

## `d6ecd534` translation groups

| Group | Upstream paths | CatBLST disposition |
| --- | --- | --- |
| Selected vector primitives | `cpu/vector/{avx512ifma,conversion,multiplication}.hpp` | Translate only operations required by BLS12-381; preserve `_mm`/`_vv`/`_cc` separation |
| Portable emulation | `cpu/vector/fallback.hpp` | Do not copy; implement equivalent operations in CatBLST's existing `cc` emulation |
| Precomputation | `cpu/precompute/precompute.hpp` | Validation/generation reference only; no production runtime computation |
| MatrixNoK parameters | `src/bls12_381_{intrns4_system,perm,perm_batch,perm_params,ring_params}.hpp` | Import exact BLS12-381 constants with provenance and independent validation |
| RNS core | `src/{bounded_element,bounded_ring,bounds,change_base,elementwise_reduce,matrix_nok_batch,ring_element,rns_reduce_perm}.hpp` | Translate required algorithms and bound types into `blst/avx/rns` |
| Conversion/inversion | `src/conversion_inversion.hpp` | Translate required conversion path; retain CatBLST inverse contracts |
| Tower | `src/{fp2,fp12}.hpp` | Port only deltas required by selected ring; preserve CatBLST field APIs |
| Curve | `src/ec.hpp` | Translate the BLS12-381 complete G1 formulas and GLV integration through CatBLST's typed candidate; preserve CatBLST point contracts and reject research-only API expansion |
| BN254 | `src/bn254_*` and BN254 branches in generic files | reject |
| Bench/tests/artifact | benchmark, test, scripts, results, examples, GPU, vendored BLST | reject as product code; use selected vectors/results only as external evidence |

## Required evidence before closing the translation row

- Exact path-by-path semantic diff against the baseline import.
- Independent validation of every selected modulus, permutation, inverse, and
  correction constant.
- Local portable-AVX Debug differential tests for each exposed operation.
- Peach native AVX-512IFMA correctness and code-generation inspection.
- One separately approved serial performance comparison against the current
  CatBLST backend.
- Formal proofs and residual assumptions recorded in `PROOF_PLAN.md`.

Confirmed defects and the evidence needed for eventual upstream contributions
are tracked separately in `VROOM_UPSTREAM_FINDINGS.md`.
