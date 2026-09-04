# BTX Prepared-Pairing Range Evidence

## Purpose

BTX verifies one ciphertext-side pairing product before decryption shares
arrive. The existing AVX implementation converts both operands to affine form
and prepares every G2 operand during each evaluation. At 256 ciphertexts, that
operation takes about 13 ms on the deployment-class Peach host. A monolithic
13 ms job also delays higher-priority work in the planned BTX executor.

This record evaluates a narrower CatBLST interface. The caller prepares the G2
operands once, assigns one immutable index to each G1 operand, and evaluates
any contiguous logical range. Each range returns a Miller-loop product before
final exponentiation. The caller multiplies the range products and applies one
final exponentiation to the aggregate.

The measured candidate is based on CatBLST commit
`d587c01df47978898f4e54ba7e180b3f65013640`. The project owner approved the
implementation and evidence batch on 2026-09-03.

## Interface and trust boundary

The C operation is
`blst_avx_prepared_g2_batch_miller_loop_product_indices`. Its inputs are an
immutable prepared-G2 batch, one native G1 slice, and one equally sized index
slice. Each index selects the prepared G2 operand paired with the corresponding
G1 point.

The operation has these properties:

- It retains no cursor, scheduler, cancellation flag, or mutable workspace.
- It allocates no memory.
- Independent workers may share the immutable prepared batch.
- Each invocation uses fixed stack arrays for at most 16 finite pairs.
- Identity operands contribute the multiplicative identity.
- Empty input returns the multiplicative identity in the Miller domain.
- Rejected null pointers or out-of-range indices leave a nonnull output
  unchanged.
- No C++ exception crosses the C boundary.
- Timing may depend on public lengths, indices, and point identities.

The caller must provide library-produced, subgroup-validated G1 points. Batch
construction receives validated native G2 points and prepares their Miller
lines. CatBLST does not repeat subgroup checks in the evaluation hot path.

Rust exposes the result as `MillerLoopResult`, not `Gt`. This distinction
prevents a pre-final Miller value from entering target-group operations. The
type exposes only identity construction, multiplication with another Miller
product, and final exponentiation.

BTX already retains the prepared G2 data used by alpha processing. That table
stores operands in reverse order. The index slice lets ciphertext verification
reuse that table with the mapping `index = bmax - 1 - position`. BTX therefore
does not need a second prepared table or a representation conversion.

## Environment and method

- Host: Peach, AMD Ryzen 9 7950X, with SMT disabled.
- Compiler: Clang and Clang++ 21.1.8.
- Build: Cargo release, `-O3 -march=native`, and Rust
  `-C target-cpu=native`.
- Backend: native AVX-512/IFMA.
- Serial affinity: CPU 8.
- Parallel affinity: CPUs 8 through 11.
- Isolation: no two benchmark processes ran concurrently.
- Sampling: Criterion used 500 ms warm-up intervals and two-second requested
  measurement intervals. Product cases used ten samples. Range cases used 20
  samples.

The benchmark keeps native points, the prepared table, index vectors, Rayon
pools, and the legacy workspace outside the timed region. The unprepared case
includes affine conversion and G2 line preparation because the current BTX
operation performs that work on every invocation. Every product case includes
exactly one final exponentiation.

The benchmark commands selected one family at a time. Serial families used
this command shape:

```text
CC=clang-21 CXX=clang++-21 CFLAGS='-O3 -march=native' \
CXXFLAGS='-O3 -march=native' RUSTFLAGS='-C target-cpu=native' \
taskset -c 8 cargo bench --manifest-path bindings/rust/Cargo.toml \
  --features avx --bench btx_prepared_pairing_product -- \
  /prepared_full/ --noplot
```

Parallel families replaced `taskset -c 8` with `taskset -c 8-11` and selected
`/parallel_/`.

## Product results

All times are Criterion center estimates in milliseconds. `Prepared-10` runs
ten-pair ranges serially, multiplies their Miller products, and applies one
final exponentiation. The parallel cases use the same ten-pair ranges and one
caller-owned Rayon pool.

| Pairs | Unprepared | Prepared full | Prepared-10 | 2 workers | 4 workers |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 32 | 1.6989 | 0.82974 | 0.85429 | 0.56833 | 0.33403 |
| 64 | 3.3230 | 1.5851 | 1.6190 | 0.92049 | 0.57764 |
| 128 | 6.5722 | 3.0960 | 3.1485 | 1.7566 | 1.0193 |
| 256 | 13.069 | 6.1186 | 6.2251 | 3.2556 | 1.7869 |

Preparing G2 lines once reduces serial evaluation time by 51.16% to 53.18%.
Ten-pair partitioning adds 2.96%, 2.14%, 1.70%, and 1.74% at the four anchors.
At 128 and 256 pairs, the serial partition cost satisfies the 2% gate.

Four workers are 2.48, 2.74, 3.04, and 3.42 times faster than one full
prepared call. Relative to the existing unprepared path, four workers are
5.09, 5.75, 6.45, and 7.31 times faster.

These parallel measurements establish kernel scalability. They do not measure
the final BTX task manager, its bounded submission window, or competition with
other BTX work.

## Range duration

The range benchmark measures one Miller product and excludes aggregation and
final exponentiation.

| Pairs per range | Estimate | 95% confidence interval |
| ---: | ---: | ---: |
| 8 | 194.57 us | 194.54--194.60 us |
| 9 | 217.01 us | 216.97--217.07 us |
| 10 | 239.44 us | 239.41--239.47 us |
| 11 | 261.93 us | 261.88--261.97 us |
| 12 | 284.42 us | 284.38--284.47 us |
| 16 | 374.07 us | 374.04--374.11 us |

Ten pairs is the largest measured grain whose center estimate is below 250 us.
The confidence interval describes the estimated mean; it is not a p95 latency
measurement. BTX integration must measure the p95 duration of scheduled jobs
before claiming the executor's 250 us latency target.

The current evidence supports ten pairs as the provisional execution grain.
The 32- and 64-pair anchors miss the 2% serial-overhead gate by 0.96 and 0.14
percentage points. The task-manager benchmark must decide whether to retain
parallel ranges for those anchors or submit one prepared monolithic job.

## Correctness and hardening evidence

The Rust tests compare the new operation with direct pairings. They cover:

- complete and empty products;
- G1 and G2 identities;
- arbitrary prepared-index permutations;
- partition grains 1, 2, 5, 10, 16, 31, and 37;
- reverse range-completion order;
- concurrent calls that share one immutable prepared batch;
- mismatched slice lengths; and
- out-of-range prepared indices.

The raw C adapter checks one direct pairing and failure-atomic rejection of an
invalid index. Peach compiled that adapter with Clang 21 under ASan and UBSan.
The native AVX corpus replay passed without a sanitizer finding.

The AVX unit and scalar/AVX conformance tests pass. The broad Cargo command has
one unrelated existing failure: rustdoc parses a C-style example in generated
`bindings.rs` as Rust. The implementation tests complete before that doctest.

## Decision

The prepared-G2 range primitive is suitable for CatBLST. It is a general
allocation-free kernel, reuses existing prepared storage, and materially
improves serial and parallel performance.

The CatBLST batch does not select BTX scheduling policy. After approval, Monad
can pin the CatBLST commit and evaluate ten-pair jobs through the real task
manager. That integration must measure job p95, end-to-end preprocessing time,
and interference with release-critical tasks.
