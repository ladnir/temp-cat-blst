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
Internal infinity tests retain exact-zero checks because library-produced
infinity has an exact zero coordinate by contract.

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

Completed locally:

- production static library compilation under portable-AVX Debug;
- all 22 focused MatrixNoK scalar-differential, curve, GT, FFT, pairing,
  production-alias, and redundant-equality tests;
- independent constant regeneration and corruption self-test;
- generated Lean snapshot reproducibility;
- focused `GtGeneratorConstant.lean` check; and
- structural production-dispatch validation.

Pending permission-gated evidence:

- relink and execute the broader C API/FFT slices in the BTx-linked main test
  executable;
- Peach Clang 21 native AVX-512IFMA integration and code-generation review;
- accept `AVX-P4-002`, currently exact-hash-bound and `in_review`, after its
  native evidence is recorded; and
- any broad suite, sanitizer, or benchmark gate.

No CUDA, BTx/BTE, benchmark, portable-AVX Release, or broad build belongs to
this packet.
