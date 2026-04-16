#pragma once

#include <span>
#include <vector>

#include "BTE/beat/types.h"

namespace bte::beat
{
    SetupResult Setup(const BeatParameters& params, const RandomBytesSource& rng);
    PreparedPublicParameters PreparePublicParameters(const PublicParameters& public_parameters);

    EncryptionResult Encrypt(
        const BeatParameters& params,
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

    BatchDecryptionPrecomputation PrecomputeBatchDecryption(
        const PublicParameters& public_parameters,
        std::span<const Ciphertext> ciphertexts);

    BatchDecryptionPrecomputation PrecomputeBatchDecryption(
        const PreparedPublicParameters& prepared_public_parameters,
        std::span<const Ciphertext> ciphertexts);

    std::vector<Bytes> FinalizeBatchDecryption(
        const BeatParameters& params,
        const PublicParameters& public_parameters,
        const BatchDecryptionPrecomputation& precomputation,
        const blst::G1Point& aggregate_predecryption_key);

    std::vector<Bytes> FinalizeBatchDecryption(
        const BeatParameters& params,
        const PreparedPublicParameters& prepared_public_parameters,
        const BatchDecryptionPrecomputation& precomputation,
        const blst::G1Point& aggregate_predecryption_key);

    std::vector<Bytes> RecoverBatchKeys(
        const BeatParameters& params,
        const PublicParameters& public_parameters,
        std::span<const Ciphertext> ciphertexts,
        const blst::G1Point& aggregate_predecryption_key);

    std::vector<Bytes> RecoverBatchKeys(
        const BeatParameters& params,
        const PreparedPublicParameters& prepared_public_parameters,
        std::span<const Ciphertext> ciphertexts,
        const blst::G1Point& aggregate_predecryption_key);
}
