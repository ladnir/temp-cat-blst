# `blst/avx/rns/matrixnok_fp.h` Review

## Scope

This record covers the complete retained BLS12-381 MatrixNoK base-field owner
through commit `d69be47`.

## Contracts

- `Element` stores one field value in both fixed eight-lane RNS bases. Its
  physical and RNS bounds are compile-time producer claims.
- `from_canonical_limbs` accepts an ordinary canonical integer in `[0,p)`;
  `from_scalar` accepts a scalar-backend Montgomery value; and
  `from_canonical_50bit_digits` accepts eight radix-`2^50` digits whose
  reconstructed ordinary integer is in `[0,p)`. Debug asserts enforce these
  narrow import preconditions without adding Release work.
- The tracked represented-multiple interval is `[-1932,2377]`. Completion adds
  `1932*p`, translating the interval to `[0,4309]`. The resulting wide
  redundant envelope fits strictly inside the exact product of the two RNS
  bases.
- One-shot and fixed-batch multiplication use the selected exact-bound
  MatrixNoK kernels. Delayed reduction remains only for formulas that combine
  several products before reduction.
- `convert_from_raw` requires distinct input/high/low vector objects and emits
  raw radix-`2^52` accumulators for the scalar conversion bridge, not
  independently reduced field values.
- Fixed batches are nonempty. Bound reattachment preserves the carrier's exact
  `rns_bounds_type`.

## Findings

- `AUD-253`: specialization retained the completed ceiling `4309` but dropped
  its derivation from VROOM's `MIN_RNS=-1932`, `MAX_RNS=2377`, and combined
  two-basis capacity check.
- `AUD-254`: generic negative-bound preparation directly negated its signed
  template lower bound without first constraining it to the reviewed range.
- `AUD-255`: an internal runtime hex parser and legacy GLV/GLS constructors
  duplicated generated compile-time constants and could silently import zero
  from malformed literals in Release.
- `AUD-256`: the raw export bridge hardcoded eight lanes and did not state its
  disjoint-object or raw-output contract at the field owner.
- `AUD-257`: bound reattachment retained a stale cv-stripping workaround, and
  the private M2 batch-expansion boundary did not state its nonempty contract.
- `AUD-264`: the comparison-only Fp2 microkernel retained two private
  base-field completion/negation hooks and privileged friendship after the
  selected generic Fp2 owner became the sole production implementation.
- `AUD-354`: the fixed-array production batch API used an `int` non-type
  parameter even though `std::array` encodes its extent as `std::size_t`, so
  GCC could not deduce the otherwise valid compile-time batch.

## Validation

- `KernelBounds.lean` proves `4309 = 2377 - (-1932)`, the completion offset
  `1932`, and the exact inequality
  `4309 * 40^2 * p^2 < M1 * M2`.
- The independent constants validator checks the same completed range against
  the generated bases, and the dispatch validator pins all repaired source
  contracts and rejects the deleted parser and constructors.
- The serial portable-AVX Debug production/test target compiles. Seven focused
  canonical import/export, base multiplication, GLV/GLS, and point-equality
  tests pass.
- All 20 accepted AVX correspondence rows and 59 mapped hashes validate after
  rebinding the Fp and scalar-multiplication owners.
- The comparison-only Fp2 hooks and friendship are absent, while selected
  one-shot, fixed-batch, delayed-reduction, conversion, and expansion callers
  continue to compile and pass their focused differentials.
- `AUD-354` aligns the production `batch_modmul` non-type parameter with
  `std::array`'s `std::size_t` extent. GCC 13 syntax-checks the exact
  portable-AVX test translation unit and MSVC compiles the focused Debug object;
  the generated operations and runtime schedule are unchanged.

## Residual risk

The compile-time bounds remain producer claims rather than runtime lane
validation. The later `AUD-259`/`AUD-263` conversion-owner review independently
bounds byte weaving through REDC and inspects MSVC carry lowering. Scalar REDC
semantics, native Clang lowering, and x86 IFMA semantics remain external. No
native code-generation run, benchmark, broad suite, CUDA, portable-AVX Release
build, or BTx/BTE audit was performed for this original Fp-owner packet.
