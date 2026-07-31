#pragma once

#include <span>
#include <vector>

#include "BTE/btx/types.h"

namespace bte::btx
{
    SetupResult Setup(const BtxParameters& params, const RandomBytesSource& rng);
    PreparedPublicParameters PreparePublicParameters(const PublicParameters& public_parameters);

    EncryptionResult Encrypt(
        const BtxParameters& params,
        const PublicParameters& public_parameters,
        const RandomBytesSource& rng);

    ThresholdRequest BuildThresholdRequest(
        const PublicParameters& public_parameters,
        std::span<const Ciphertext> ciphertexts);

    PartialDecryptionShare CreatePartialDecryptionShare(
        const ServerSecretShare& secret_share,
        const ThresholdRequest& request);

    blst::G1Point CombinePartialDecryptionShares(
        std::span<const PartialDecryptionShare> shares,
        const PublicParameters& public_parameters);

    bool VerifyAggregatePredecryptionKey(
        const ThresholdRequest& request,
        const blst::G1Point& aggregate_predecryption_key,
        const PublicParameters& public_parameters);

    bool VerifyAggregatePredecryptionKey(
        const ThresholdRequest& request,
        const blst::G1Point& aggregate_predecryption_key,
        const PublicParameters& public_parameters,
        AggregatePredecryptionKeyVerificationWorkspace& workspace);

    std::vector<blst::GT> ComputeCrossTerms(
        const PublicParameters& public_parameters,
        std::span<const Ciphertext> ciphertexts,
        CrossTermAlgorithm algorithm = CrossTermAlgorithm::MiddleProductFft);

    std::vector<blst::GT> ComputeCrossTerms(
        const PreparedPublicParameters& prepared_public_parameters,
        std::span<const Ciphertext> ciphertexts,
        CrossTermAlgorithm algorithm = CrossTermAlgorithm::MiddleProductFft);

    BatchDecryptionPrecomputation PrecomputeBatchDecryption(
        const PublicParameters& public_parameters,
        std::span<const Ciphertext> ciphertexts,
        CrossTermAlgorithm algorithm = CrossTermAlgorithm::MiddleProductFft);

    BatchDecryptionPrecomputation PrecomputeBatchDecryption(
        const PreparedPublicParameters& prepared_public_parameters,
        std::span<const Ciphertext> ciphertexts,
        CrossTermAlgorithm algorithm = CrossTermAlgorithm::MiddleProductFft);

    std::vector<Bytes> FinalizeBatchDecryption(
        const BtxParameters& params,
        const PublicParameters& public_parameters,
        const BatchDecryptionPrecomputation& precomputation,
        const blst::G1Point& aggregate_predecryption_key);

    std::vector<Bytes> FinalizeBatchDecryption(
        const BtxParameters& params,
        const PreparedPublicParameters& prepared_public_parameters,
        const BatchDecryptionPrecomputation& precomputation,
        const blst::G1Point& aggregate_predecryption_key);

    std::vector<Bytes> RecoverBatchKeys(
        const BtxParameters& params,
        const PublicParameters& public_parameters,
        std::span<const Ciphertext> ciphertexts,
        const blst::G1Point& aggregate_predecryption_key,
        CrossTermAlgorithm algorithm = CrossTermAlgorithm::MiddleProductFft);

    std::vector<Bytes> RecoverBatchKeys(
        const BtxParameters& params,
        const PreparedPublicParameters& prepared_public_parameters,
        std::span<const Ciphertext> ciphertexts,
        const blst::G1Point& aggregate_predecryption_key,
        CrossTermAlgorithm algorithm = CrossTermAlgorithm::MiddleProductFft);
}
