# BTX Scalar Prepared-Pairing Range Evidence

## Purpose

The scalar BTX backend needs the same independently schedulable
ciphertext-side pairing product as the native AVX backend. Repeating G2 affine
conversion and Miller-line preparation in every evaluation wastes work, while
one monolithic scalar pairing product prevents other BTX jobs from using the
application's worker pool.

This batch adds one immutable scalar prepared-G2 owner and one stateless
indexed range operation. It does not add a scheduler or a second protocol
algorithm. Rust may partition one product into ranges, evaluate those ranges
concurrently, multiply the pre-final products as a tree, and apply final
exponentiation once.

The review candidate is based on CatBLST commit `295c17c4`. It remains
uncommitted pending project-owner approval.

## Interface and trust boundary

`blst_scalar_prepared_g2_batch_new` accepts a nonempty sequence of
library-produced, subgroup-validated scalar G2 projective points. Construction
allocates one immutable 68-line table and one identity byte per operand. G2
identity is valid and is recorded without reading an uninitialized line table.
Allocation or construction failure returns null, and no C++ exception crosses
the C boundary.

`blst_scalar_prepared_g2_batch_miller_loop_product_indices` accepts native
scalar G1 projective points and equally many prepared indices. It validates all
indices before arithmetic. Each call then:

- removes G1 or G2 identity contributions;
- collects at most 16 finite G1 operands in fixed stack storage;
- converts the finite G1 operands to affine form with one batched inversion;
- evaluates the matching prepared G2 lines with the scalar 16-way Miller
  schedule; and
- multiplies chunk products into one pre-final Miller result.

The operation allocates no memory, retains no cursor or mutable workspace, and
may execute concurrently against one immutable prepared batch. Empty input
returns the Miller multiplicative identity. Invalid pointers or indices leave
a nonnull output unchanged. Timing may depend on the public count, indices,
and identity locations. The caller establishes subgroup membership once at the
typed or decoding boundary; the hot operation does not repeat it.

Rust exposes the owner as `scalar::PreparedG2Batch` and the result as
`MillerLoopResult`. Shape and index errors are typed and are checked before
FFI. RAII releases the native table. `Send + Sync` is justified by immutable
native storage and const-only evaluation.

## Environment and method

- Host: Peach, AMD Ryzen 9 7950X, SMT disabled.
- Compiler: Clang and Clang++ 21.1.8.
- Native flags: release, `-O3 -march=native`; Rust
  `-C target-cpu=native`.
- Affinity: CPUs 8 through 11 for serial and four-worker comparisons.
- Isolation: no benchmark process ran concurrently with another benchmark.
- Inputs: nonidentity scalar multiples `17*G1` and `29*G2`, retained in native
  projective representation outside the timed region.
- Sampling: Criterion used 500 ms warm-up, a two-second requested measurement
  interval, and ten samples for product cases.

The unprepared comparison converts both G1 and G2 operands to affine form,
evaluates the existing multi-Miller loop, and final-exponentiates. Prepared
cases keep the line table outside the timed region. Every complete product case
includes exactly one final exponentiation.

## Product results

Times are Criterion center estimates in milliseconds. `Serial ranges` uses
16-pair calls and multiplies their Miller products serially. `Parallel` uses
four worker threads, eight-pair ranges at 32 pairs, and 16-pair ranges at all
larger anchors.

| Pairs | Unprepared | Prepared full | Serial ranges | Parallel | Prepared speedup |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 32 | 4.9732 | 2.9850 | 2.9792 | 1.0144 | 2.94x |
| 64 | 9.6614 | 5.6901 | 5.6809 | 1.6604 | 3.43x |
| 128 | 19.048 | 11.099 | 11.086 | 3.0240 | 3.67x |
| 256 | 37.803 | 21.907 | 21.891 | 5.7417 | 3.82x |

Preparation removes 40.0% to 42.1% of repeated serial evaluation time.
Sixteen-pair serial partitioning measures 0.07% to 0.19% faster than the full
prepared call; that difference is noise, but it establishes near-zero range
overhead at the four anchors.

At 32 pairs, two 16-pair jobs use only two workers and take 1.6971 ms. Four
eight-pair jobs take 1.0144 ms. For 64, 128, and 256 pairs, 16-pair jobs are
faster than eight-pair jobs: 1.6604 versus 1.7361 ms, 3.0240 versus 3.1738 ms,
and 5.7417 versus 6.0762 ms. The initial scalar policy should therefore derive
grain from the logical count and available worker slots, while capping a range
at the native 16-pair schedule.

## Range duration and scheduling tradeoff

One two-pair range measures 232.12 us. One 16-pair range measures 1.3493 ms.
The latter amortizes affine inversion and preserves the native Miller batch;
forcing every range below 250 us would discard most of that throughput.

The 250 us target is therefore not a valid hard bound for this scalar phase.
The BTX executor must still expose parallel work and priority at range
boundaries. It should select the smallest grain that keeps the available
workers occupied without breaking the 16-pair batching advantage, then record
the resulting nonpreemptible duration. Release-critical work can overtake
unsubmitted pairing ranges; a submitted scalar range runs to completion.

## Correctness and hardening evidence

The safe Rust tests cover complete and empty products, arbitrary prepared
index permutations, partitions of 1, 2, 5, 10, 16, 31, and 37 pairs, reverse
completion order, concurrent sharing, G1 and G2 identities, length mismatch,
and out-of-range indices. The raw C test covers direct equality, split-product
composition, identity, null handling, and failure-atomic rejection.

On Peach:

- all 73 Rust library tests pass in Clang 21 release/native mode;
- the full Clang 21 release/native C++ suite passes 580 tests and skips four
  unconfigured external-vector tests;
- the focused scalar prepared-batch test passes under ASan and UBSan; and
- evaluation succeeds with global C++ allocation deliberately disabled.

The complete sanitizer suite stops on the pre-existing null-output probe in
`blst_scalar_from_hexascii`; the focused new path is clean. This batch does not
expand the broad API fuzzer. Deterministic C and Rust tests exercise successful
construction/evaluation and failure-atomic invalid-index rejection directly.

## Decision

The scalar prepared-G2 range primitive is suitable for CatBLST. It is native,
allocation-free during evaluation, immutable across workers, failure-atomic,
and materially faster than repeated preparation. The range contract supports
BTX's required parallel path without placing scheduling policy in C++.

This evidence does not complete BTX scheduling. Monad must still combine range
products with a parallel reduction, select grain using its actual pool width,
and include this operation in the end-to-end linear-work completion matrix.
