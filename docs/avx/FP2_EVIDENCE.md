# Selected MatrixNoK Fp2 Evidence

## Scope

- Selected implementation audited through: `13fe98a`
- Upstream algorithm commit: `d6ecd534c54a41702ea80ae75b52250b6d51b372`
- Compared VROOM head: `26bd9aee2168f167583cf096add355547a7b56d4`
- C++ owner: `blst/avx/fields/fp2.h`
- Base-field owner: `blst/avx/rns/matrixnok_fp.h`
- Formal owner: `proof/avx/CatBLST/AVX/MatrixNoK/Fp2.lean`

Production uses `FP2Ring<FpRingCandidate>`. Values represent `x + u*y`, where
`u^2 = -1`. `prep_left` prepares `(x,y,-y)` once so delayed multiplication
forms `ac-bd` and `ad+bc` without signed accumulators. Squaring forms
`a^2-b^2` and `2ab`. Reduction batches all real coordinates before all
imaginary coordinates.

The owner is allocation-free and fixed-schedule with respect to field data.
It stores the empty fixed base ring by value with `no_unique_address`, so the
selected owner has the one-byte empty-object minimum and cannot contain a
dangling ring reference. Default initialization remains overwrite-before-read;
value initialization zeroes the selected trivial carriers. `is_zero` and
`equal` compare exact RNS storage and are not semantic field observations when
redundant representatives are possible.

## Correctness and bounds

Accepted correspondence row `AVX-P3-001` proves:

- prepared negation represents `-y` in both RNS bases;
- multiplication and squaring implement the Fp2 formulas;
- `mul_accumulate`, `mul_sub_accumulate`, base-field scaling, the nonresidue,
  and `3b` formulas are correct;
- the strided batch layout reconstructs each original Fp2 value;
- each two-product coordinate fits the tight `8m²` accumulator bound;
- completion, M1 reduction, batched R1 change-base, M2 reduction, and R2
  correction return the exact selected representative; and
- the completion-boundary witness remains inside the proved range.

`AvxMatrixNoKPrimitives.SelectedFp2MultiplyAndSquareMatchScalar` compares the
selected implementation with scalar `blst::Fp2::mul` and `blst::Fp2::sqr` for
zero, one, `p-1`, dense values, all 64 pairings of eight Fp2 values, every
square, assignment back to either source object, a 64-product dependency
chain, a 16-square dependency chain, and the `<2m_i` output bounds.
`SelectedFp2CompletionBoundaryMatchesScalar` retains the weak R2-margin
witness. The fixed-context scalar-accumulator oracle also passes.

## Selected optimized code generation

`test/audit_avx_matrixnok_fp2_codegen.cpp` is an intentionally small
compile-only probe. It avoids instantiating the complete curve and pairing
stack when only selected Fp2 code generation is under review.

MSVC 19.50 compiled it with C++23, `/O2`, `/Ob3`, `/DNDEBUG`, and
`/arch:AVX512`:

| Operation | Instructions | Frame | Calls | Conditional branches |
| --- | ---: | ---: | ---: | ---: |
| multiply | 371 | 1,264 bytes | 0 | 5 |
| square | 328 | 1,024 bytes | 0 | 5 |

All five branches in each body are fixed-control matrix/permutation loops;
none depends on field data. No Fp2, `_vv`, reduction, or owner wrapper remains
out of line. The ownership and contract repairs change no arithmetic
expression or operation schedule.

## Removed comparison microkernel

The former `matrixnok::Fp2RingCandidate` was an isolated comparison owner. It
was never dispatched and was outside the formal correspondence even though it
used the same field formula. Its private base-ring hooks, tests, and probe
duplicated selected functionality after integration, so `AUD-264` removes it
at `13fe98a`. The selected scalar differential and code-generation probe now
own the retained evidence.

Historical Clang 21 measurements of that deleted microkernel were:

| Operation | Prior CatBLST | Deleted microkernel | Pinned VROOM |
| --- | ---: | ---: | ---: |
| multiply | 41.4 ns | 33.6 ns | 34.7 ns |
| square | 41.1 ns | 33.6 ns | 34.1 ns |

Those figures explain why the experimental implementation was originally
retained. They are provenance only and must not be presented as measurements
of the selected delayed implementation. No benchmark ran for this audit.

## Exact correspondence inputs

- `fp2.h` SHA-256:
  `a065c9ccc4e691d6322b1220c9cf56ba07ddbbb126364fd3cb59ad8e18ea47f9`
- `matrixnok_fp.h` SHA-256:
  `ef9524d0895b1e072e15a8a8e9950723ddd4d76372d70e687f62397bbec82060`

The hash update is the `AUD-354` `std::size_t` fixed-batch interface repair;
the selected Fp2 arithmetic and generated code contract are unchanged.

The canonical exact-hash ledger is `proof/avx/CORRESPONDENCE.csv`. Higher Fp12,
GT, curve, and pairing claims remain owned by their separate accepted rows.

## Residual risk

x86 IFMA semantics and compiler lowering remain the shared P0 boundary.
Milestone sanitizer, fuzzing, and statistical leakage evidence are outside
this file review. No remote run, broad suite, CUDA, BTx/BTE audit, or
portable-AVX Release build occurred in the `13fe98a` packet. Commit `65206c9`
adds only the trivial default constructor used for value-initialized prepared
line storage; selected Fp2 arithmetic and its accepted theorems are unchanged.
