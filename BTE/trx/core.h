#pragma once

#include <optional>
#include <span>
#include <vector>

#include "BTE/trx/types.h"

namespace bte::trx
{
    Crs GenerateCRS(const TrxParameters& params, const RandomBytesSource& rng);

    ValidatorPublicShare DeriveValidatorPublicShare(
        const ValidatorSecretShare& secret_share,
        const blst::G2Point& h_tau);

    PublicKeyMaterial DerivePublicKeyMaterial(
        std::span<const ValidatorPublicShare> public_shares,
        std::span<const std::uint32_t> subset_indices);

    EncryptionKey ToEncryptionKey(const PublicKeyMaterial& public_key_material);

    EncryptedCiphertext Encrypt(
        const TrxParameters& params,
        const EncryptionKey& encryption_key,
        std::span<const std::uint8_t> associated_data,
        const RandomBytesSource& rng);

    EncryptedCiphertext Encrypt(
        const TrxParameters& params,
        const EncryptionKey& encryption_key,
        const RandomBytesSource& rng);

    std::optional<VerifiedCiphertext> VerifyCiphertext(
        const TrxParameters& params,
        const Ciphertext& ciphertext,
        std::span<const std::uint8_t> associated_data);

    std::optional<VerifiedCiphertext> VerifyCiphertext(
        const TrxParameters& params,
        const Ciphertext& ciphertext);

    BatchDigest BuildBatchDigest(
        const TrxParameters& params,
        std::span<const CiphertextEnvelope> ciphertexts,
        std::size_t context_id,
        const Crs& crs);

    BatchDigest BuildBatchDigest(
        const TrxParameters& params,
        std::span<const Ciphertext> ciphertexts,
        std::size_t context_id,
        const Crs& crs);

    std::vector<blst::G1Point> ComputeEvalProofs(
        const BatchDigest& digest,
        const Crs& crs);

    bool VerifyEvaluationProof(
        const BatchDigest& digest,
        std::size_t proof_index,
        const blst::G1Point& proof,
        const Crs& crs);

    PartialDecryptionShare CreatePartialDecryption(
        const TrxParameters& params,
        const ValidatorSecretShare& secret_share,
        const BatchDigest& digest,
        const PublicKeyMaterial& public_key_material);

    bool VerifyPartialDecryption(
        const TrxParameters& params,
        const PartialDecryptionShare& share,
        const ValidatorPublicShare& public_share,
        const BatchDigest& digest,
        const PublicKeyMaterial& public_key_material);

    blst::G1Point CombineShares(
        std::span<const PartialDecryptionShare> shares,
        const BatchDigest& digest,
        const PublicKeyMaterial& public_key_material);

    std::optional<Bytes> DecapsulateOne(
        const TrxParameters& params,
        const BatchDigest& digest,
        std::size_t ciphertext_index,
        const blst::G1Point& proof,
        const blst::G1Point& sigma,
        const Crs& crs,
        const PublicKeyMaterial& public_key_material);
}
