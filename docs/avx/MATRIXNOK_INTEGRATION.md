# MatrixNoK Production Integration

## Selected call graph

Production AVX dispatch now selects the fixed MatrixNoK construction through
semantic aliases rather than using `Candidate` names at call sites:

- `matrixnok::FpRing` owns Fp;
- `matrixnok::G2Fp2Ring` owns Fp2 and G2;
- `matrixnok::Fp12Ring` owns Fp12 and GT;
- `matrixnok::G1Curve`, `G2Curve`, `G1GLV`, and `G2GLS` own curve arithmetic;
- `matrixnok::MillerEngine` and `PairingInverter` own pairing; and
- G1/GT FFT and middle-product adapters use the same selected field owners.

The exact `Candidate` spellings remain the source/theorem correspondence
owners. The production aliases are compile-time type aliases: they introduce
no runtime dispatch, conversion, allocation, or wrapper call. The incumbent
`BoundedRing` implementation remains in-tree temporarily for historical
differential/reference work, but it is unreachable from the production point,
GT, pairing, and GT-FFT owners. `validate_matrixnok_dispatch.py` fails closed
if those owners regain an incumbent type or stale `52`-specific adapter name.

## Integration repairs

MatrixNoK standard elements have redundant representatives. A field zero may
therefore be encoded by several nonzero raw lane vectors. The inherited point
equality path treated raw-lane zero as semantic field zero, producing false
inequality after valid arithmetic. Projective equality now keeps its native
cross-products and canonicalizes only the final two differences at the
observation boundary. Affine equality canonicalizes its coordinates through
the same scalar bridge. This adds work only to explicit equality/validation;
curve arithmetic, GLV/GLS, MSM, FFT, Miller, and pairing kernels are unchanged.
Internal arithmetic owners may retain exact-zero checks only where their local
representation invariant proves exact zero. Public infinity/identity
observations and conversion boundaries use canonical field-zero observation.
The native integration run showed why this distinction is required: valid
low-order multiplication can produce a redundant nonzero RNS lift of zero.
Secret projective normalization derives its zero bit from the scalar export
already required for inversion, so identity signing remains constant-time and
does not add a second conversion.

The old GT-generator literal belonged to the incumbent RNS bases. Production
now uses a fixed MatrixNoK literal generated from CatBLST's canonical scalar GT
serialization. The independent Python validator reconstructs all 192 lanes;
`GtGeneratorConstant.lean` proves their exact standard encoding and lane
contracts. The assertion that the twelve canonical Fp coefficients are
CatBLST's chosen GT generator remains an explicit reviewed scalar-constant
boundary and a scalar/AVX differential test.

The opaque pairing context continues to begin object lifetimes with placement
construction. MatrixNoK Fp2 storage is trivially default constructible while
value-initialization still produces zero, so the C boundary performs no hidden
initialization work and cannot throw during context construction.

## Evidence status

Completed locally and on Peach:

- production static library compilation under portable-AVX Debug;
- all 22 portable-Debug MatrixNoK scalar-differential tests;
- a fresh Clang 21.1.8 native Release production build at `6ce0d44`;
- all 34 native focused MatrixNoK, C API contract, and native-point contract
  tests, including redundant-zero infinity and identity-signing regressions;
- independent constant regeneration and corruption self-test;
- generated Lean snapshot reproducibility;
- focused `GtGeneratorConstant.lean` check; and
- structural production-dispatch validation;
- all three native MatrixNoK audit CTests; and
- optimized native code-generation inspection with no out-of-line `_vv` or
  conversion symbols and no stack-frame change from the preceding commit.

Pending permission-gated evidence:

- relink and execute the broader C API/FFT slices in the BTx-linked main test
  executable;
- any broad suite, sanitizer, or benchmark gate.

`AVX-P4-002` is accepted at the exact hashes recorded in the correspondence
ledger. No CUDA, BTx/BTE, benchmark, portable-AVX Release, or broad build
belongs to this packet.

## Lean integration-gap audit

Production selection does not turn the existing proof into a blanket proof of
the entire AVX API. The exact current boundary is:

- proved and source-bound: IFMA/portable primitives, MatrixNoK conversion and
  reduction, Fp/Fp2/Fp12 arithmetic, GT scalar schedules, final exponentiation,
  and the retained direct/prepared Miller loop;
- newly closed: the production GT literal's 192-lane standard encoding;
- structurally and natively bound: the zero-overhead production aliases in
  accepted `AVX-P4-002`;
- formally proved over the selected projective groups: public G1/G2 complete
  point addition, including infinity, doubling, inverse, and generic inputs;
- formally proved over the conventional cyclic G1/G2 domains: public masking,
  modulo-`r` preparation, concrete generator/eigenbasis facts, and all four
  selected CT/VT GLV/GLS scalar-multiplication schedules;
- partially closed formal evidence: the exact relation-based public membership
  predicates now satisfy
  `standard cyclic subgroup ⊆ accepted relation ⊆ r-torsion` without an
  ambient-cardinality premise; identifying the accepted rational eigenspace
  with the standard generator span remains `AUD-231`;
- explicit non-Lean observation boundary: MatrixNoK canonical export through
  scalar REDC/byte weaving for point equality and validation;
- outside this arithmetic proof: generic FFT orchestration, map/hash-to-curve
  protocol composition, C ABI state machines, compiler/ISA correctness, and
  microarchitectural constant-time behavior; and
- retained earlier residual: the concrete Fp12 inversion provider used by
  final exponentiation is differential evidence rather than a Lean inversion
  theorem.

These are evidence gaps rather than confirmed C++ defects. `AUD-227` is closed
for selected public curve arithmetic by accepted row `AVX-P4-003`; proof
commit `24de7b4` discharges relation soundness and cyclic-generator inclusion,
while `AUD-231` keeps the final eigenspace/span claim visible so future release
notes cannot imply broader coverage than the correspondence ledger provides.

Native object inspection also found function-local dynamic initialization
guards around several derived AVX owners and public literals. `AUD-229` is
validated at `fe09e26`: selected owners are constant-initialized, public C accessors
return the immutable native literals directly, and an independently regenerated
curve-constant artifact replaces conversion during initialization. The change
does not alter the arithmetic kernels or the established generated constants.
Local MSVC object inspection finds no dynamic-init symbols. Peach Clang 21.1.8
Release passes 35 focused tests and all five audits; guard/atexit relocations
fall to zero, the selected arithmetic probe is byte-identical, and hot owner
stack frames are unchanged. `AUD-230` canonicalizes exact proof hashes to LF,
closing the Windows-CRLF/Peach-LF audit-infrastructure defect exposed by this
run. No benchmark was warranted because native code generation was unchanged.
