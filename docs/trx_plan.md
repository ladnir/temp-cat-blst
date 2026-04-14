# TrX Core Implementation Plan

## Scope

This document captures the planned implementation of the plain TrX batched threshold encryption core in this repository, adapted to the current `blst` codebase and project structure.

This plan is for:

- C++ implementation
- built on the existing `blst` pairing stack
- no consensus integration
- no Tess / witness-encryption path
- external DKG for threshold shares
- trusted CRS generation
- parameterized batch size and context count

This plan is not for:

- consensus or mempool hooks
- epoch policy
- non-interactive key generation
- full payload AEAD design beyond encapsulated session-key recovery

## Key Decisions

### One-time signature choice

We will use a one-time BLS signature built from `blst`, not Ed25519.

We will use the `min-pk` layout:

- one-time verification key in `G1`
- one-time signature in `G2`
- `vk_sig` serialized in compressed `G1` form, 48 bytes
- `sig` serialized in compressed `G2` form, 96 bytes

This keeps the verification key small while staying entirely inside the existing `blst` stack.

### Group conventions

We will use the same group placement as the TrX construction:

- KZG commitments and evaluation proofs in `G1`
- threshold public key and ciphertext base terms in `G2`
- recovered encapsulation secret in `GT`

### Project layout

Instead of forcing TrX into the current main `blst` target immediately, we will create a separate `BTE` project area and build its crypto core on top of `blst`.

Recommended layout:

- `BTE/`
- `BTE/CMakeLists.txt`
- `BTE/trx/`
- `BTE/building_blocks/`
- `test/`

Recommended files:

`BTE/trx/`

- `params.h`
- `types.h`
- `crs.h`
- `ciphertext.h`
- `digest.h`
- `threshold.h`
- `core.h`

`BTE/building_blocks/`

- `polynomial.h`
- `kzg.h`
- `lagrange.h`
- `serialization.h`
- `hkdf_sha256.h`
- `bls_ots.h`

`test/`

- `Trx_core_tests.cpp`
- `Trx_negative_tests.cpp`
- `Trx_vectors_tests.cpp`

## Repo Fit

## What already exists

The current repo already contains the main primitives needed for TrX:

- scalar field arithmetic in `Fr`
- `G1`, `G2`, and `GT`
- compressed serialization and deserialization
- subgroup and infinity checks
- pairings, including Miller loop and final exponentiation
- MSM support
- FFT and polynomial-oriented infrastructure
- SHA-256
- internal HKDF code
- BLS signing and verification

Important existing surfaces:

- pairing and GT operations: `blst/scalar/curve/GT.h`, `blst/scalar/curve/GT.cpp`
- point serialization and subgroup checks: `blst/scalar/curve/point.h`
- C++ wrapper surface: `bindings/blst.hpp`
- MSM helpers: `blst/curve/multi_scalar.h`
- FFT and coefficient contexts: `blst/fft/ctFft.h`, `blst/fft/context.h`
- SHA-256 and internal HKDF code: `blst/sign/sha256.h`, `blst/sign/keygen.cpp`

## What needs adaptation

The repo does not currently expose a reusable public HKDF helper for arbitrary callers. The HKDF logic exists, but today it is embedded inside BLS keygen code.

That means TrX does not need a new cryptographic primitive for KDF, but it does need one of:

- refactoring existing HKDF into a shared helper, or
- adding a BTE-local HKDF-SHA256 helper that reuses the same SHA backend

The existing `middle_product` path is also not the direct implementation vehicle for TrX KZG commitments and proofs. It is specialized around pairing-style public convolutions. We should reuse ideas and lower-level helpers from that stack, but implement TrX KZG operations directly over `Fr` polynomials and `G1` CRS powers.

## Cryptographic Objects

## Parameters

Define a parameter object containing:

- `max_batch_size`
- `num_contexts`
- domain separators for `H1`, `HF`, KDF, and serialization
- serialization version
- point encoding policy

## CRS

For each context `dc`, sample:

- global `tau in Fr*`
- context randomizer `kappa_dc in Fr*`

Store:

- `contexts[dc].g1_powers[j] = g^(kappa_dc * tau^j)` for `0 <= j <= max_batch_size`
- `h_tau = h^tau`

The CRS must include degree `0` through `max_batch_size`.

## Threshold key material

TrX consumes externally generated threshold shares:

- secret share `sk_i`
- validator index `i`
- public share `pk_i = h^(sk_i)`
- tau-share public element `pk_tau_i = (h^tau)^(sk_i)`

From any threshold subset, reconstruct:

- aggregate `pk = h^sk`
- aggregate `pk_tau = h^(tau * sk)`

using Lagrange interpolation in the exponent.

## Ciphertext Design

## One-time BLS signing key

Generate a one-time BLS `min-pk` keypair:

- secret signing key `sk_sig`
- verification key `vk_sig in G1`

Serialize `vk_sig` in compressed `G1` form.

## Tag

Associated data is carried explicitly through encryption, verification, and
batch construction.

Define:

- `tg = HF(vk_sig_bytes, associated_data)`

where `HF : bytes -> Fr`.

## Encapsulation

Sample:

- random session key bytes `K_session`
- random nonzero `alpha in Fr`

Compute:

- `ct1 = (pk_tau * pk^(-tg))^alpha` in `G2`
- `ct2 = h^alpha` in `G2`
- `S = e(H1(pk), pk)^alpha` in `GT`
- `mask = KDF_GT(S)`
- `ct3 = K_session xor mask`

Sign the canonical ciphertext header with the one-time BLS key.

Ciphertext fields:

- `ct1`
- `ct2`
- `ct3`
- `vk_sig`
- `sig`

## Verification

Before a ciphertext enters a batch:

1. deserialize and subgroup-check `ct1`
2. deserialize and subgroup-check `ct2`
3. deserialize and subgroup-check `vk_sig`
4. reject unexpected infinity values
5. verify the one-time BLS signature on the canonical ciphertext header
6. recompute `tg = HF(vk_sig_bytes, associated_data)`

Invalid ciphertexts are excluded before batch polynomial formation.

## Batch Digest

Given valid ciphertexts for one context:

1. compute tags `tg_i`
2. reject duplicate tags in v1
3. construct the vanishing polynomial

`f(X) = product_i (X - tg_i)`

4. commit with the selected context CRS

`com = g^(kappa_dc * f(tau))`

Implementation:

- represent `f` by coefficients in `Fr`
- compute `com` as a `G1` MSM against the context power vector

Digest output should cache:

- valid batch order
- tags
- `f` coefficients
- `com`
- batch size
- context id

## Partial Decryption

For validator share `sk_i`:

- `base = H1(pk) * com^(-1)` in `G1`
- `pd_i = base^(sk_i)` in `G1`

Store the share together with the validator index.

## Share Verification

Verify:

`e(H1(pk) * com^(-1), pk_i) == e(pd_i, h)`

using the exact same batch digest and commitment object that will be used in combination.

## Share Combination

For threshold subset `S`, compute Lagrange coefficients at `0` from the explicit validator indices:

`lambda_i = product_{j in S, j != i} (-j)/(i-j)`

Then combine:

`sigma = product_{i in S} pd_i^(lambda_i)`

This yields:

`sigma = (H1(pk) * com^(-1))^sk`

Duplicate indices must be rejected.

## Evaluation Proofs

For each ciphertext tag `tg_i`:

- compute `q_i(X) = f(X)/(X - tg_i)`
- compute `pi_i = g^(kappa_dc * q_i(tau))`

Implemented strategy:

- build the full vanishing polynomial once
- derive all exclusion quotients with a divide-and-conquer product-tree path
- do one `G1` MSM per proof

This removes the synthetic-division bottleneck while keeping the simpler
one-commitment-per-proof structure.

## Decapsulation

For ciphertext `i`:

- `S_i = e(pi_i, ct1_i) * e(sigma, ct2_i)` in `GT`
- `mask_i = KDF_GT(S_i)`
- `K_session_i = ct3_i xor mask_i`

## Building Blocks

## Hashes

### `H1`

Map serialized aggregate `pk` to `G1` using proper hash-to-curve with domain separation.

Since this repo already supports hash-to-curve on `G1`, this is a natural fit.

### `HF`

Map serialized one-time `vk_sig` bytes to `Fr` with a TrX-specific DST.

This should reduce a hash output modulo `r` in a deterministic, documented way.

### `KDF_GT`

Serialize `GT` canonically, then feed the bytes into HKDF-SHA256 or an equivalent SHA-256-based KDF.

Implementation note:

- the primitive exists in this repo already
- the shared API does not yet exist
- BTE should either factor it out of `blst/sign/keygen.cpp` or provide a local helper

## Canonical serialization

All signed and hashed structures must use one canonical encoding:

- compressed points
- fixed byte order
- explicit version byte
- length prefixes where needed

Never sign raw in-memory structs.

## How BTE should use existing blst code

Use directly:

- `Fr`, `G1Point`, `G2Point`, `GT`
- point serialization and subgroup checks
- `GT::miller_loop`, `GT::final_mult`
- existing MSM primitives
- FFT roots and coefficient contexts where helpful

Do not treat `middle_product` as the TrX KZG engine. Reuse its lessons, not its full API.

Use the internal public C++ APIs for the one-time BLS signature path.

## Recommended API Surface

- `GenerateCRS(max_batch_size, num_contexts)`
- `DerivePublicKeyMaterial(dkg_public_shares, dkg_tau_shares, subset)`
- `Encrypt(ek, associated_data, rng)`
- `VerifyCiphertext(ciphertext, associated_data)`
- `BuildBatchDigest(ciphertext_envelopes, context_id, crs)`
- `ComputeEvalProofs(batch_digest, crs)`
- `CreatePartialDecryption(sk_share, batch_digest, pk)`
- `VerifyPartialDecryption(pd_i, pk_i, batch_digest, pk)`
- `CombineShares(valid_shares, subset_indices, batch_digest, pk)`
- `DecapsulateOne(ciphertext_i, pi_i, sigma, pk)`

## Testing Plan

## Unit tests

- field and point serialization round-trips
- subgroup and infinity checks
- pairing bilinearity sanity checks
- MSM versus naive multiplication
- KZG commit and proof correctness
- threshold interpolation in the exponent

## Protocol tests

- end-to-end happy path
- threshold boundary: `t` succeeds, `t-1` fails
- batch-size sweep
- context isolation
- deterministic vectors

## Negative tests

- malformed `ct1`
- malformed `ct2`
- malformed `vk_sig`
- bad one-time BLS signature
- altered tag
- wrong commitment
- wrong validator index
- shuffled proof
- wrong `pk_tau`
- duplicate share
- duplicate tag
- wrong context

## Implementation Sequence

### Stage 1

- create BTE structure
- define TrX types and serialization
- implement CRS container
- implement DKG boundary types
- implement one-time BLS signature helper
- implement HKDF-SHA256 helper

### Stage 2

- implement encryption and verification
- implement batch digest construction
- implement partial decryption, share verification, and combination
- implement divide-and-conquer quotient generation plus proof commitments
- implement decapsulation

### Stage 3

- add deterministic vectors
- add negative tests
- benchmark commitments and proofs
- optimize proof generation

## Current conclusion

With the updated constraints, the TrX plan fits this repo well.

The main repo-specific adjustments are:

- use one-time BLS signatures in `min-pk` mode instead of Ed25519
- treat HKDF as an existing primitive that needs a reusable helper surface
- implement TrX as a separate BTE project area that builds on `blst`
- reuse `blst` algebra, MSM, FFT, serialization, and pairing code directly

The pairing implementation is present and usable. The main engineering work is integration, API shaping, canonical encoding, and the new TrX/KZG/threshold logic itself.

## Backend Selection

BTE stays written against the top-level `blst::*` C++ API.

That means backend choice is static and happens at configure/build time through `BLST_DEFAULT_AVX`:

- `BLST_DEFAULT_AVX=OFF` selects the scalar backend
- `BLST_DEFAULT_AVX=ON` selects the AVX backend

The same BTE source builds in either mode because the top-level aliases in `blst/curve/point.h`, `blst/curve/GT.h`, and `blst/curve/map_to_G.h` resolve to the configured backend.

For convenience, the repository now provides these Windows debug presets:

- `x64-Debug-Scalar`
- `x64-Debug-Avx`

Example commands:

- `cmake --preset x64-Debug-Scalar`
- `cmake --build out/build/x64-Debug-Scalar`
- `ctest --test-dir out/build/x64-Debug-Scalar --output-on-failure`

- `cmake --preset x64-Debug-Avx`
- `cmake --build out/build/x64-Debug-Avx`
- `ctest --test-dir out/build/x64-Debug-Avx --output-on-failure`

This is a compile-time choice, not a runtime switch. If both backends are needed, build BTE twice with different presets.

## Implemented V1 Notes

The current repository state now includes a working v1 BTE/TrX core implementation on top of `blst`.

Implementation choices made during execution:

- the BTE implementation now stays on the C++ side of the repo surface: the TrX core and its helpers use the repo's internal `blst` C++ field, point, hash-to-curve, and pairing types directly
- HKDF is implemented as a BTE-local SHA-256 helper in `BTE/building_blocks/hkdf_sha256.*`
- duplicate tags are rejected during batch construction
- ciphertext verification returns an optional verified object and filters invalid ciphertexts before digest construction
- randomness enters BTE as a byte-filling callback, and CRS/session scalars are derived from full-width random seed bytes through a C++ hash-to-scalar helper rather than `uint64_t` truncation
- evaluation-proof verification is implemented with the pairing-based zero-opening relation
- the scalar-to-byte path for point and GT multiplication must convert `Fr` exactly once; the implemented code now uses `Fr::to_scaler` directly, which fixed a double-conversion bug that otherwise broke interpolation and decapsulation

Test interpretation choices made during execution:

- context-specific commitments, proofs, and combined partial decryptions differ across contexts
- a fully self-consistent proof-plus-sigma pair from another context can still decrypt the same ciphertext batch because the context randomizer cancels algebraically
- the enforced isolation property in tests is therefore that mixing artifacts across contexts fails, not that a self-consistent alternate-context artifact set necessarily fails

## Paper Alignment Actions

The following items were identified by comparing the BTE implementation against
the Section 4 / Section 6 TrX construction in the paper and have now been
implemented in this pass.

### Correctness hardening

- associated data is now carried explicitly through `Encrypt`,
  `VerifyCiphertext`, `BuildBatchDigest`, and the signed ciphertext header so
  tags follow `HF(vkSig, ad)`
- the checked decryption path now enforces the paper's validity gates:
  - combine rejects duplicate indices and fewer-than-threshold shares
  - decapsulation rejects invalid evaluation proofs
- negative coverage was extended for:
  - malformed `ct2`
  - wrong validator index
  - shuffled proof
  - wrong commitment
  - wrong `pk_tau`

### Performance work

- `ComputeEvalProofs` now uses a divide-and-conquer all-quotients path based on
  the product tree
- `ComputeEvalProofs` now uses an FK-style recursive proof engine that derives
  child proof CRS slices by group/scalar convolution over the product tree
- the BTE/TrX path has been validated under both static backend presets:
  - `x64-Debug-Scalar`
  - `x64-Debug-Avx`

### Deferred work after this pass

- broader caching of digest-side trees and common interpolation subsets
- dedicated microbenchmarks for TrX subroutines beyond the stage-level harness

## Current Status

The current BTE/TrX implementation now includes:

- associated-data-aware tags and signed ciphertext headers
- checked share combination with threshold enforcement
- checked decapsulation with proof verification
- fast `O(n log^2 n)` Lagrange interpolation
- faster all-quotients evaluation-proof preparation
- cached `H1(pk)` reuse in key material and cached evaluation proofs on the digest
- deterministic regression vectors for the single-ciphertext happy path
- expanded negative coverage for the paper-alignment cases above

Backend validation completed successfully for the dedicated BTE test target:

- default backend: `out/build/x64-Debug/test/blst_tests.exe --gtest_filter=Trx.*`
- scalar backend: `out/build/x64-Debug-Scalar/test/bte_trx_tests.exe --gtest_filter=Trx.*`
- AVX backend: `out/build/x64-Debug-Avx/test/bte_trx_tests.exe --gtest_filter=Trx.*`

Full scalar-suite validation also passes:

- `cmake --build out/build/x64-Debug-Scalar --target blst_tests --config Debug`
- `ctest --test-dir out/build/x64-Debug-Scalar --output-on-failure --timeout 1200`

## Benchmark Findings

The repository now includes a dedicated BTE/TrX benchmark harness in
`benchmark/Trx_bench.cpp`, wired through the `bte_trx_benchmarks` target.

The benchmark is organized around the major TrX stages:

- `GenerateCRS`
- `DerivePublicKeyMaterial`
- `Encrypt`
- `VerifyCiphertext`
- `BuildBatchDigest`
- `ComputeEvalProofs`
- `CreatePartialDecryption`
- `VerifyPartialDecryption`
- `CombineShares`
- `DecapsulateOne`
- `EndToEndBatch`

The `ComputeEvalProofs` benchmark explicitly clears the digest proof cache per
iteration so it measures the uncached FK proof-engine cost rather than the
cached fast-return path.

Current measured scalar results on the local machine show:

- `ComputeEvalProofs` is still the dominant bottleneck, but much smaller than
  before the FK-style proof-engine pass
- `EndToEndBatch` is dominated by proof generation, with per-ciphertext
  decapsulation next
- `BuildBatchDigest` is materially cheaper than proof generation
- share creation and share combination are comparatively small costs

Representative scalar measurements at batch size `128`:

- `BuildBatchDigest/128`: about `0.20s` CPU
- `ComputeEvalProofs/128`: about `1.23s` CPU
- `DecapsulateOne/128`: about `2.60ms` per ciphertext
- `EndToEndBatch/128`: about `1.86s` CPU

VTune hotspot runs were collected for:

- `BuildBatchDigest/128`
- `ComputeEvalProofs/128`
- `EndToEndBatch/128`

The profiles show that scalar Montgomery field arithmetic dominates these stages,
especially `mul_mont_3841`. This indicates that the implementation is currently
limited by the expected algebraic work rather than by control-flow overhead.

## Updated Performance Conclusion

The current implementation has already removed the old synthetic-division
bottleneck and now uses:

- fast vanishing-polynomial construction
- fast `O(n log^2 n)` interpolation
- divide-and-conquer all-quotients generation
- an FK-style recursive proof engine that derives proof commitments through
  product-tree recursion plus specialized middle-window `G1 x Fr -> G1`
  convolution

This significantly reduces proof-generation time in practice, and changes the
performance picture from “proof generation overwhelmingly dominates” to “proof
generation is still the largest stage, but by a much smaller margin.”

In practical terms:

- the main missing single-threaded performance work is now further caching and
  tuning on top of the FK-style proof engine
- parallelism is intentionally deferred for now
- the next single-threaded optimization areas are caching and MSM tuning, not a
  major algorithmic redesign

## Remaining Performance Gaps

The dominant remaining cost is still inside the FK proof engine's child-CRS
derivation. In the current implementation, each recursive FK node derives child
CRS slices by doing a full group/scalar convolution and then extracting the
needed coefficient window.

The main gaps identified in this path are:

- child-CRS derivation was still structurally a middle-window convolution, but
  it was not using a dedicated middle-product helper
- the generic FFT path was paying for work outside the needed output window
- the FK path was not fully aligned with the repo's truncated FFT utilities
  that already exist in `blst/fft/ctFft.h`

## Next Concrete Optimization

The next concrete optimization is to treat FK child-CRS derivation explicitly as
a middle-product problem:

- input coefficients in `G1`
- sibling polynomial coefficients in `Fr`
- output coefficients in `G1`
- only the middle window of the convolution is needed

Concretely, the helper should:

- reverse the sibling polynomial once
- run truncated forward FFTs using the true nonzero prefix lengths
- do the pointwise `G1 x Fr -> G1` multiplication in the evaluation domain
- use truncated inverse FFT to recover only the required middle window

This is a specialized `G1 x Fr -> G1` middle-window convolution helper. It is
in the same algorithmic family as the repo's middle-product machinery, but it is
not the same public API because the existing `middle_product(...)` entrypoints
are specialized around `G1 x G2 -> GT` pairings and multipair code.

This optimization has now been implemented in the BTE KZG layer. FK child-CRS
derivation now uses a dedicated `G1 x Fr -> G1` middle-window helper built on:

- truncated forward FFT with true nonzero prefix lengths
- pointwise `G1 x Fr -> G1` multiplication in the evaluation domain
- truncated inverse FFT over only the required coefficient window
- per-call caches for FFT roots and truncated FFT plans

Follow-up scalar benchmarking and VTune on this path show that the remaining
cost is now concentrated in the low-level arithmetic used by the FFT and
group/scalar multiplication kernels, not in a missing higher-level algorithm.
The dominant hotspots remain:

- `mul_mont_3841`
- `mul_3841`
- `sub_mod_384`
- `add_mod_384`

This means the next tuning targets are below the TrX algorithm level:

- improving the FFT kernel shape used by the `G1 x Fr -> G1` helper
- exploiting more twiddle special cases inside the group/scalar FFT path
- reducing the cost of repeated `G1` scalar multiplications inside the
  evaluation-domain multiplication and inverse-range recovery

A generic twiddle fast-path experiment for `0`, `1`, and `-1` was also tried in
the shared FFT templates, but it did not produce a reliable scalar win on the
current FK proof path, so it was not kept. Any future twiddle optimization
should likely be done as a more specialized kernel for the `G1 x Fr -> G1`
middle-window path rather than as a generic template-layer branch.

An attempt was also made to push all FK child-CRS scratch buffers into a
preallocated per-depth workspace in the current one-shot `ComputeEvalProofs`
API. That change was not kept either. The implementation was correct, but it
did not produce a stable scalar win in the uncached proof benchmark, and it
complicated the code path substantially. If we want true "no dynamic
allocations in the hot loop" semantics, the cleaner design is a separate
reusable FK workspace object that can be prepared once and reused across proof
computations, rather than forcing that preparation into the default one-shot
API.
