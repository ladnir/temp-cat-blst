#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "BTE/transparent_bte/types.h"

namespace bte::transparent_bte
{
    /**
     * Build deterministic Shamir material from caller-supplied coefficients.
     * This is a benchmark/test adapter, not a DKG. `sharing_polynomial[0]` is
     * the secret and the vector length is the threshold.
     */
    SetupResult SetupFromSharingPolynomial(
        std::span<const blst::Fr> sharing_polynomial,
        std::size_t num_servers,
        const blst::G1Point& independent_generator);

    /** Encrypt the supplied GT message with caller-supplied nonzero randomness. */
    Ciphertext Encrypt(
        const PublicParameters& public_parameters,
        const blst::Fr& randomness,
        const blst::GT& message);

    /** Compute H(c1) exactly as used by the Cauchy transform. */
    blst::Fr CiphertextPoint(const Ciphertext& ciphertext);

    /**
     * Produce this server's single share for a nonempty admitted batch.
     * The caller verifies ciphertext proofs and pairwise-distinct H(c1) values
     * once at admission; this algebraic experiment does not implement the
     * paper's uninstantiated NIZK.
     */
    PartialDecryptionShare PartialDecrypt(
        const SecretShare& secret_share,
        std::span<const Ciphertext> ciphertexts);

    /**
     * Verify one share against a nonempty admitted batch. Returns false for an
     * invalid share index or pairing equation. Ciphertext admission is the
     * caller's responsibility and is intentionally not repeated here.
     */
    bool VerifyPartialDecryptionShare(
        const PublicParameters& public_parameters,
        const PartialDecryptionShare& share,
        std::span<const Ciphertext> ciphertexts);

    blst::G1Point CombinePartialDecryptionShares(
        const PublicParameters& public_parameters,
        std::span<const PartialDecryptionShare> shares);

    /**
     * Direct O(B^2) implementation used as a correctness oracle.
     * Throws when the batch is empty or two H(c1) values collide.
     */
    BatchPrecomputation PrecomputeBatchDecryptionNaive(
        std::span<const Ciphertext> ciphertexts);

    /**
     * Direct O(B^2) Cauchy transform using one native public-scalar G1 MSM and
     * one native public-scalar G2 MSM per output row. The AVX implementation
     * operates on native projective points and native Fr values throughout.
     */
    BatchPrecomputation PrecomputeBatchDecryptionMsm(
        std::span<const Ciphertext> ciphertexts);

    /** Benchmark-only scalar distribution control for the native row-MSM path. */
    enum class MsmWeightMode
    {
        Cauchy,
        Unit,
    };

    /** Phase timings produced by the benchmark-only row-MSM entry point. */
    struct MsmPrecomputeTimings
    {
        std::uint64_t point_hash_ns = 0;
        std::uint64_t output_shell_ns = 0;
        std::uint64_t cauchy_weights_ns = 0;
        std::uint64_t buffer_setup_ns = 0;
        std::uint64_t g1_msm_ns = 0;
        std::uint64_t g2_msm_ns = 0;
        std::uint64_t affine_ns = 0;
        std::uint64_t pairing_ns = 0;
        std::uint64_t total_ns = 0;
    };

    /**
     * Benchmark-only form of `PrecomputeBatchDecryptionMsm`. `Unit` still
     * constructs the real Cauchy weights, then replaces them with one before
     * entering the otherwise identical MSM and pairing path. Its output is not
     * a protocol result. Scalar builds support only `Cauchy`.
     */
    BatchPrecomputation BenchmarkPrecomputeBatchDecryptionMsm(
        std::span<const Ciphertext> ciphertexts,
        MsmWeightMode weight_mode,
        MsmPrecomputeTimings* timings);

    /**
     * Direct O(B^2) prepared-pairing implementation. It reuses every Cauchy
     * pair's scaled G1 inputs across both affected output slots, precomputes
     * each G2 line table once, and performs one final exponentiation per slot.
     * Throws when the batch is empty or two H(c1) values collide.
     */
    BatchPrecomputation PrecomputeBatchDecryptionQuadraticPairing(
        std::span<const Ciphertext> ciphertexts);

    /**
     * Lemma 1 subproduct-tree/remainder-tree implementation.
     * Throws when the batch is empty or two H(c1) values collide. All scalar
     * multiplications in this variable-time path operate on public values.
     */
    BatchPrecomputation PrecomputeBatchDecryptionFast(
        std::span<const Ciphertext> ciphertexts);

    /**
     * Finish Eq. (1) from a combined share. Returns one GT plaintext per slot.
     * The precomputation and ciphertext spans must describe the same batch.
     */
    std::vector<blst::GT> FinalizeBatchDecryption(
        const BatchPrecomputation& precomputation,
        const blst::G1Point& combined_share,
        std::span<const Ciphertext> ciphertexts);
}
