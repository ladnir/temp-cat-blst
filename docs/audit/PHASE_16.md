# Phase 16: AVX Curves and Protocols

## Scope

Phase 16 performs the file-by-file hand audit of the native and portable AVX
curve, pairing, map-to-group, and multi-scalar implementation. Scalar code is
used as an independent differential oracle where appropriate. Native AVX-512
code generation and runtime evidence is required before performance-sensitive
files close; portable-AVX Debug supplies fast local correctness feedback.

BTx/BTE remains frozen, CUDA remains experimental, and neither is audited by
this phase. Portable-AVX Release, broad suites, sanitizers, fuzzing, remote
machines, and benchmarks remain explicit milestone or permission gates.

## AVX point review

The owning `blst/avx/curve/point.h` review follows the completed selected
facade review. It traces the AVX homogeneous representation, scalar boundary
conversions, fixed-window exact multiplication, GLV/GLS multiplication,
signing and key derivation, basic point operations, validation, encoding, and
visible helper functions through native C++, the selected C API, and focused
tests.

Initial findings are:

- `AUD-169`: signing and public-key derivation select AVX types but execute
  scalar point multiplication and redundant normalizations.
- `AUD-170`: equality and zero observation use scalar conversions and
  inversions; the dormant cross-product helper lacks infinity parity.
- `AUD-171`: `to_pk2` declarations lost the inherited pre-decorated-output
  contract and in AVX C++ state the wrong size and format.
- `AUD-172`: `times_minus_z` requires distinct input and output objects but did
  not state or assert that invariant.
- `AUD-173`: AVX selection removes the inherited native C++ `sign2` method and
  forces the C boundary to duplicate normalization and output logic.

The approved remediation keeps multiplication and coordinate arithmetic
AVX-native, retains a fixed-schedule scalar base-field inversion only as the
normalization boundary until an independently justified native inversion
exists, removes observation conversions, restores exact output documentation,
and adds zero-Release-overhead Debug invariants.

## Validation policy

The edit loop uses only the existing selected portable-AVX Debug target and
structural documentation checks. Focused tests cover G1/G2 scalar agreement,
canonical and all-zero infinity representatives, projective scaling, signing,
public-key derivation, optional outputs, pre-decorated finalization, exact
aliasing where promised, and rejected aliasing where prohibited.

Native AVX-512/IFMA runtime and optimized code-generation inspection remain
required before `AUD-169` and the file review are validated. No benchmark runs
without explicit approval.
