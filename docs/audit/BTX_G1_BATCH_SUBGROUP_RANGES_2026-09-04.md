# BTX G1 Batch-Subgroup Range Evidence — 2026-09-04

## Decision

CatBLST now exposes a stateless range kernel for its existing randomized G1
batch-subgroup check. The misuse-resistant one-shot validator remains the
default validation boundary. BTX may use the range kernel only after canonical
decoding has established that every native point is on-curve and only through
a driver that checks every challenge row exactly once.

CatBLST owns challenge construction. A fresh 32-byte seed expands through the
existing domain-separated SHA-256 stream into 81 independent ternary rows. For
the selected curve, a fixed invalid batch passes all rows with probability at
most `(1/3)^81 < 2^-128` in the SHA-256 PRF model. Parallel execution does not
change this argument: ranges partition the same fixed challenge and the caller
accepts only the conjunction of every range result.

## Interface and trust boundary

The C interface has three operations:

- a checked byte-size query for the complete row-major challenge;
- failure-atomic challenge initialization from caller-provided entropy; and
- an allocation-free half-open row-range check.

The range call does not repeat curve checks. Every point must be nonnull,
on-curve, immutable, and resident in the native AVX representation until all
ranges complete. The ordered pointer array and challenge buffer are immutable
shared inputs. Each call has no output buffer, retained state, cursor, callback,
allocation, or scheduler. CatBLST validates the range, complete challenge
capacity, nonnull pointers, and coefficient alphabet before returning success.

Rust remains responsible for owning those buffers, preserving their lifetime,
partitioning all 81 rows, applying the stage barrier, and combining results.
An aggregate rejection never identifies a bad point; the Rust driver retains
the exact per-ciphertext fallback that filters invalid candidates.

## Correctness and hardening evidence

Focused native AVX tests cover:

- exact size, empty input, and overflow;
- undersized buffers and entropy failure without output publication;
- exception containment at the entropy callback boundary;
- deterministic challenge expansion and coefficient range;
- disjoint row partitions and invalid/reversed ranges;
- malformed coefficients confined to their selected range; and
- rejection of a known on-curve order-11 G1 point across irregular partitions.

The complete 17-test AVX C-API contract suite passes under Clang 21 Release.
The two new tests also pass in a combined native-plus-portable Clang 21 build
under AddressSanitizer and UndefinedBehaviorSanitizer.

## Peach performance evidence

Host: AMD Ryzen 9 7950X. Compiler: Clang 21.1.8. Build: Release,
`-O3 -DNDEBUG -march=native`, native AVX-512, assembly and ADX enabled. Each
entry is the median of five sequential repetitions. No other benchmark process
ran concurrently.

The ranged benchmark includes challenge generation, creation and joining of
fresh `std::thread` workers on every iteration, and the conjunction of all
range results. It is therefore conservative relative to BTX's persistent
Rayon pool.

| Points | One-shot | Ranged 1 | Ranged 2 | Ranged 4 | Ranged 8 | Ranged 16 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 64 | 2.09 ms | 2.14 ms | 1.14 ms | 0.653 ms | 0.467 ms | 0.503 ms |
| 256 | 3.78 ms | 3.78 ms | 2.02 ms | 1.15 ms | 0.778 ms | 0.715 ms |
| 1,024 | 10.5 ms | 10.3 ms | 5.57 ms | 3.17 ms | 2.00 ms | 1.54 ms |

At 1,024 points, the one-worker ranged path is 1.9% faster than the one-shot
path and therefore has no material serial overhead. Sixteen fresh workers
reduce wall time by 85.3%, a 6.8x speedup. The 64-point result saturates at
eight workers, so the Rust scheduler should retain a size-aware sequential or
low-worker policy rather than forcing maximum parallelism at every anchor.

## Reproduction

Configure CatBLST on Peach with Clang 21, native AVX-512, tests and benchmarks,
then run:

```text
taskset -c 0-15 \
  out/benchmark/avx_g1_batch_subgroup_benchmarks \
  --benchmark_repetitions=5 \
  --benchmark_report_aggregates_only=true
```

The benchmark target is `avx_g1_batch_subgroup_benchmarks`. The production BTX
integration must be measured separately because Rayon task grain, the decode
barrier, exact fallback, and proof verification determine end-to-end admission
latency.
