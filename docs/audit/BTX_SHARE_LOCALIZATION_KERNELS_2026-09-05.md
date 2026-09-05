# BTX share-localization kernel evidence, 2026-09-05

## Purpose

BTX first checks a combined decryption share. A failed check causes the Rust
implementation to verify subsets with fresh random weights. For a batch of
256 ciphertexts, each visited subset constructs 257 public-scalar G2 sums and
evaluates 258 Miller-loop terms. This record describes the CatBLST kernels
used to reduce that cost.

The change does not define the localization algorithm. Rust still chooses the
subsets, samples the weights, evaluates the fault tree, and decides which
shares survive. CatBLST receives validated subgroup points and public reduced
scalars for one existing equation.

## Interfaces

`blst_avx_p2s_mult_fr_vartime_scratch_sizeof(max_points)` returns the reusable
scratch capacity for every input length through `max_points`.

`blst_avx_p2s_mult_fr_vartime` computes

```text
out = scalars[0] * points[0] + ... + scalars[n - 1] * points[n - 1].
```

The points are native AVX G2 projective values. They must already be on the
curve and in the prime-order subgroup. Each scalar is a public reduced
`blst_fr` value in native Montgomery form. The function is variable-time and
must not receive secret scalars.

`blst_avx_miller_loop_n` computes one unfinalized product for parallel arrays
of native affine G1 and G2 points. Rust composes that result with the remaining
terms and performs one final exponentiation for the complete fault equation.

The Rust binding exposes `G2MsmWorkspace` and `PairingWorkspace`. Construction
allocates aligned storage. Evaluation reuses that storage and allocates
nothing.

## G2 MSM algorithm

For at most 32 points, the MSM applies the existing GLS decomposition to each
point and public scalar. A width-four wNAF table represents each of the four
GLS components. One 64-bit doubling schedule is shared across all components
and all points in the selected batch.

The compiled batch sizes are 2, 3, 4, 6, 8, 12, and 16. Larger inputs through
32 are decomposed into those batches. Sizes 3, 6, and 12 match the fault tree
for a 12-share candidate set. Inputs above 32 use the existing G2 Pippenger
kernel.

The caller owns the temporary point array, wNAF tables, digit arrays, and
canonical scalar staging area. The scratch query returns the maximum required
size across the complete supported prefix. This rule matters at the 32-to-33
algorithm crossover because Pippenger can require less scratch than the
small-batch GLS kernel.

## Contract

The C operation has these preconditions:

- `out` is non-null.
- A nonempty input has non-null point and scalar arrays.
- Every point is a library-produced G2 subgroup point.
- Every scalar is public.
- Required scratch is 64-byte aligned.
- The used output, input, and scratch ranges are pairwise disjoint when the
  input contains more than one point.

Empty input writes the G2 identity. A one-point input needs no scratch and may
alias its output with its point. Invalid shape, capacity, alignment, or overlap
returns false before publishing output. The function does not allocate or
throw across the C boundary.

Concurrent calls may share immutable point and scalar inputs. Each call must
have a distinct output and distinct scratch storage. CatBLST retains no cursor,
task, priority, cancellation, or thread-pool state.

## Security scope

The optimization shares arithmetic schedules and reuses caller-owned storage.
It does not reuse verifier weights or Miller results across fault equations.
It does not change the equation, subset tree, survivor rule, or threshold.

The variable-time MSM is appropriate because fault-localization weights and
verification keys are public to the local verifier. Secret-scalar local-share
generation continues to use the constant-time G1 MSM interface.

These facts establish implementation equivalence to the existing equations.
They do not constitute a new security proof for the localization protocol.

## Performance evidence

Peach ran Clang 21 release binaries with `-O3 -march=native` and one pinned
core for the kernel benchmark. Each row is the median of 21 samples. No other
benchmark process ran concurrently.

| G2 terms | Separate scalar products | Shared-doubling MSM | MSM / separate |
|---:|---:|---:|---:|
| 2 | 40.187 us | 26.597 us | 0.662 |
| 3 | 59.697 us | 32.958 us | 0.552 |
| 4 | 80.072 us | 39.761 us | 0.497 |
| 6 | 119.178 us | 51.817 us | 0.435 |
| 8 | 158.589 us | 64.456 us | 0.406 |
| 12 | 239.704 us | 90.418 us | 0.377 |
| 16 | 318.842 us | 115.120 us | 0.361 |

The complete Monad benchmark uses 256 ciphertexts, 16 DKG indices, threshold
11, and 21 samples. At four workers, one bad share among 12 improves from
139.195 ms to 92.687 ms. One bad share among 16 improves from 167.782 ms to
105.419 ms. Five bad shares among 16 improve from 308.582 ms to 211.832 ms.
The valid-share median changes from 9.457 ms to 9.463 ms.

The final profile attributes approximately 25% of samples to batched Miller
line evaluation and 34% to the new shared-doubling G2 MSM variants. Wrapper,
allocation, and scheduler overhead do not appear among the sampled hot
symbols. Further material improvement requires pairing or G2 arithmetic work,
or a reviewed change to the localization algorithm.

## Validation

The C contract test compares independent scalar products with the public MSM
at counts 2, 3, 4, 6, 8, 12, 16, 17, 31, 32, and 33. It also covers empty
input, one-point aliasing, insufficient scratch, misalignment, overlap,
overflow, and failure atomicity.

The Rust differential test uses distinct points and scalars at the same
dispatch boundaries. A separate test covers length mismatch, empty input,
insufficient capacity, and capacity overflow. The capacity test detected the
32-to-33 sizing defect during this audit; the final scratch query is monotone
across that crossover.

The focused C tests pass under Clang 21 ASan and UBSan. The CatBLST Rust AVX
library and integration tests pass: 80 passed and one benchmark was ignored.
The two scalar/AVX conformance tests also pass. The complete AVX-default C++
run passes 358 tests and skips four vector tests whose file was not configured;
one scalar-default prepared-pairing test must run in its documented scalar
configuration.

Monad passes 173 AVX-enabled library tests and all seven shared conformance
vector tests against the final CatBLST source.
