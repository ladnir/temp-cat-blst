#include "bindings/blst_trx.h"

#include <cstring>
#include <memory>
#include <span>
#include <stdexcept>
#include <vector>

#include "BTE/trx/core.h"
#include "blst/scalar/curve/point.h"
#if BLST_DEFAULT_AVX
#include "blst/avx/curve/map_to_G.h"
#endif

namespace
{
    using BoundaryG1 = blst::scalar::G1Point;
    using BoundaryG2 = blst::scalar::G2Point;

    static_assert(sizeof(blst_p1) == sizeof(BoundaryG1));
    static_assert(alignof(blst_p1) == alignof(BoundaryG1));
    static_assert(sizeof(blst_p2) == sizeof(BoundaryG2));
    static_assert(alignof(blst_p2) == alignof(BoundaryG2));

    bte::trx::TrxParameters to_cpp_params(const blst_trx_params& params)
    {
        bte::trx::TrxParameters out{};
        out.max_batch_size = params.max_batch_size;
        out.num_contexts = params.num_contexts;
        out.session_key_size = params.session_key_size;
        out.serialization_version = params.serialization_version;
        if (params.h1_dst != nullptr) {
            out.h1_dst = params.h1_dst;
        }
        if (params.hf_dst != nullptr) {
            out.hf_dst = params.hf_dst;
        }
        if (params.kdf_salt != nullptr) {
            out.kdf_salt = params.kdf_salt;
        }
        if (params.kdf_info != nullptr) {
            out.kdf_info = params.kdf_info;
        }
        if (params.ots_sign_dst != nullptr) {
            out.ots_sign_dst = params.ots_sign_dst;
        }
        return out;
    }

    bte::trx::RandomBytesSource to_rng(
        blst_fill_random_fn fill_random,
        void* user_data)
    {
        return [fill_random, user_data](std::span<std::uint8_t> out) {
            if (fill_random == nullptr || !fill_random(user_data, out.data(), out.size())) {
                throw std::runtime_error("TRX random callback failed");
            }
        };
    }

    blst::G1Point to_public_g1(const blst_p1& point)
    {
        const auto& scalar_point = reinterpret_cast<const BoundaryG1&>(point);
#if BLST_DEFAULT_AVX
        return blst::avx::avx_map_to_G_detail::from_scalar(scalar_point);
#else
        return scalar_point;
#endif
    }

    blst::G2Point to_public_g2(const blst_p2& point)
    {
        const auto& scalar_point = reinterpret_cast<const BoundaryG2&>(point);
#if BLST_DEFAULT_AVX
        return blst::avx::avx_map_to_G_detail::from_scalar(scalar_point);
#else
        return scalar_point;
#endif
    }

    blst_p1 to_boundary_g1(const blst::G1Point& point)
    {
        blst_p1 out{};
#if BLST_DEFAULT_AVX
        const BoundaryG1 scalar_point = blst::avx::avx_map_to_G_detail::to_scalar(point);
        std::memcpy(&out, &scalar_point, sizeof(out));
#else
        std::memcpy(&out, &point, sizeof(out));
#endif
        return out;
    }

    blst_p2 to_boundary_g2(const blst::G2Point& point)
    {
        blst_p2 out{};
#if BLST_DEFAULT_AVX
        const BoundaryG2 scalar_point = blst::avx::avx_map_to_G_detail::to_scalar(point);
        std::memcpy(&out, &scalar_point, sizeof(out));
#else
        std::memcpy(&out, &point, sizeof(out));
#endif
        return out;
    }

    blst::Fr to_fr(const blst_scalar& scalar)
    {
        blst::Fr out{};
        blst::Fr::from_scalar(out, scalar.b);
        return out;
    }

    bte::trx::ValidatorSecretShare to_cpp_secret_share(
        const blst_trx_validator_secret_share& secret_share)
    {
        return bte::trx::ValidatorSecretShare{ secret_share.index, to_fr(secret_share.sk_share) };
    }

    bte::trx::ValidatorPublicShare to_cpp_public_share(
        const blst_trx_validator_public_share& public_share)
    {
        return bte::trx::ValidatorPublicShare{
            public_share.index,
            to_public_g2(public_share.pk_share),
            to_public_g2(public_share.pk_tau_share),
        };
    }

    bte::trx::PartialDecryptionShare to_cpp_partial_share(
        const blst_trx_partial_decryption_share& share)
    {
        return bte::trx::PartialDecryptionShare{ share.index, to_public_g1(share.value) };
    }
}

struct blst_trx_crs
{
    bte::trx::Crs value;
};

struct blst_trx_public_key_material
{
    bte::trx::PublicKeyMaterial value;
};

struct blst_trx_encryption_key
{
    bte::trx::EncryptionKey value;
};

struct blst_trx_encrypted_ciphertext
{
    bte::trx::EncryptedCiphertext value;
};

struct blst_trx_verified_ciphertext
{
    bte::trx::VerifiedCiphertext value;
};

struct blst_trx_batch_digest
{
    bte::trx::BatchDigest value;
};

extern "C" {

blst_trx_crs* blst_trx_generate_crs(
    const blst_trx_params* params,
    blst_fill_random_fn fill_random,
    void* user_data)
{
    if (params == nullptr || fill_random == nullptr) {
        return nullptr;
    }

    try {
        auto out = std::make_unique<blst_trx_crs>();
        out->value = bte::trx::GenerateCRS(to_cpp_params(*params), to_rng(fill_random, user_data));
        return out.release();
    } catch (...) {
        return nullptr;
    }
}

void blst_trx_crs_free(blst_trx_crs* crs)
{
    delete crs;
}

bool blst_trx_derive_validator_public_share(
    const blst_trx_validator_secret_share* secret_share,
    const blst_trx_crs* crs,
    blst_trx_validator_public_share* out_public_share)
{
    if (secret_share == nullptr || crs == nullptr || out_public_share == nullptr) {
        return false;
    }

    try {
        const bte::trx::ValidatorPublicShare public_share =
            bte::trx::DeriveValidatorPublicShare(to_cpp_secret_share(*secret_share), crs->value.h_tau);
        out_public_share->index = public_share.index;
        out_public_share->pk_share = to_boundary_g2(public_share.pk_share);
        out_public_share->pk_tau_share = to_boundary_g2(public_share.pk_tau_share);
        return true;
    } catch (...) {
        return false;
    }
}

blst_trx_public_key_material* blst_trx_derive_public_key_material(
    const blst_trx_validator_public_share* public_shares,
    const size_t public_share_count,
    const uint32_t* subset_indices,
    const size_t subset_count)
{
    if (public_shares == nullptr || public_share_count == 0 ||
        subset_indices == nullptr || subset_count == 0) {
        return nullptr;
    }

    try {
        std::vector<bte::trx::ValidatorPublicShare> cpp_public_shares;
        cpp_public_shares.reserve(public_share_count);
        for (std::size_t i = 0; i < public_share_count; ++i) {
            cpp_public_shares.push_back(to_cpp_public_share(public_shares[i]));
        }

        auto out = std::make_unique<blst_trx_public_key_material>();
        out->value = bte::trx::DerivePublicKeyMaterial(
            std::span<const bte::trx::ValidatorPublicShare>(cpp_public_shares.data(), cpp_public_shares.size()),
            std::span<const std::uint32_t>(subset_indices, subset_count));
        return out.release();
    } catch (...) {
        return nullptr;
    }
}

void blst_trx_public_key_material_free(
    blst_trx_public_key_material* public_key_material)
{
    delete public_key_material;
}

blst_trx_encryption_key* blst_trx_to_encryption_key(
    const blst_trx_public_key_material* public_key_material)
{
    if (public_key_material == nullptr) {
        return nullptr;
    }

    try {
        auto out = std::make_unique<blst_trx_encryption_key>();
        out->value = bte::trx::ToEncryptionKey(public_key_material->value);
        return out.release();
    } catch (...) {
        return nullptr;
    }
}

void blst_trx_encryption_key_free(
    blst_trx_encryption_key* encryption_key)
{
    delete encryption_key;
}

blst_trx_encrypted_ciphertext* blst_trx_encrypt(
    const blst_trx_params* params,
    const blst_trx_encryption_key* encryption_key,
    const byte* associated_data,
    const size_t associated_data_len,
    blst_fill_random_fn fill_random,
    void* user_data)
{
    if (params == nullptr || encryption_key == nullptr || fill_random == nullptr) {
        return nullptr;
    }

    try {
        auto out = std::make_unique<blst_trx_encrypted_ciphertext>();
        out->value = bte::trx::Encrypt(
            to_cpp_params(*params),
            encryption_key->value,
            std::span<const std::uint8_t>(associated_data, associated_data_len),
            to_rng(fill_random, user_data));
        return out.release();
    } catch (...) {
        return nullptr;
    }
}

size_t blst_trx_encrypted_ciphertext_session_key_size(
    const blst_trx_encrypted_ciphertext* encrypted_ciphertext)
{
    if (encrypted_ciphertext == nullptr) {
        return 0;
    }
    return encrypted_ciphertext->value.session_key.size();
}

bool blst_trx_encrypted_ciphertext_copy_session_key(
    const blst_trx_encrypted_ciphertext* encrypted_ciphertext,
    byte* out_session_key,
    const size_t out_session_key_len)
{
    if (encrypted_ciphertext == nullptr || out_session_key == nullptr) {
        return false;
    }
    if (out_session_key_len != encrypted_ciphertext->value.session_key.size()) {
        return false;
    }
    std::memcpy(
        out_session_key,
        encrypted_ciphertext->value.session_key.data(),
        encrypted_ciphertext->value.session_key.size());
    return true;
}

void blst_trx_encrypted_ciphertext_free(
    blst_trx_encrypted_ciphertext* encrypted_ciphertext)
{
    delete encrypted_ciphertext;
}

blst_trx_verified_ciphertext* blst_trx_verify_ciphertext(
    const blst_trx_params* params,
    const blst_trx_encrypted_ciphertext* encrypted_ciphertext)
{
    if (params == nullptr || encrypted_ciphertext == nullptr) {
        return nullptr;
    }

    try {
        const auto verified = bte::trx::VerifyCiphertext(
            to_cpp_params(*params),
            encrypted_ciphertext->value.ciphertext,
            std::span<const std::uint8_t>(
                encrypted_ciphertext->value.associated_data.data(),
                encrypted_ciphertext->value.associated_data.size()));
        if (!verified.has_value()) {
            return nullptr;
        }
        auto out = std::make_unique<blst_trx_verified_ciphertext>();
        out->value = *verified;
        return out.release();
    } catch (...) {
        return nullptr;
    }
}

void blst_trx_verified_ciphertext_free(
    blst_trx_verified_ciphertext* verified_ciphertext)
{
    delete verified_ciphertext;
}

blst_trx_batch_digest* blst_trx_build_batch_digest(
    const blst_trx_params* params,
    const blst_trx_verified_ciphertext* const* verified_ciphertexts,
    const size_t verified_ciphertext_count,
    const size_t context_id,
    const blst_trx_crs* crs)
{
    if (params == nullptr || verified_ciphertexts == nullptr ||
        verified_ciphertext_count == 0 || crs == nullptr) {
        return nullptr;
    }

    try {
        std::vector<bte::trx::CiphertextEnvelope> ciphertexts;
        ciphertexts.reserve(verified_ciphertext_count);
        for (std::size_t i = 0; i < verified_ciphertext_count; ++i) {
            if (verified_ciphertexts[i] == nullptr) {
                return nullptr;
            }
            ciphertexts.push_back({
                verified_ciphertexts[i]->value.ciphertext,
                verified_ciphertexts[i]->value.associated_data,
            });
        }

        auto out = std::make_unique<blst_trx_batch_digest>();
        out->value = bte::trx::BuildBatchDigest(
            to_cpp_params(*params),
            std::span<const bte::trx::CiphertextEnvelope>(ciphertexts.data(), ciphertexts.size()),
            context_id,
            crs->value);
        return out.release();
    } catch (...) {
        return nullptr;
    }
}

size_t blst_trx_batch_digest_valid_ciphertext_count(
    const blst_trx_batch_digest* batch_digest)
{
    if (batch_digest == nullptr) {
        return 0;
    }
    return batch_digest->value.valid_ciphertexts.size();
}

bool blst_trx_compute_eval_proofs(
    const blst_trx_batch_digest* batch_digest,
    const blst_trx_crs* crs,
    blst_p1* out_proofs,
    const size_t out_proof_count)
{
    if (batch_digest == nullptr || crs == nullptr || out_proofs == nullptr) {
        return false;
    }

    try {
        const std::vector<blst::G1Point> proofs = bte::trx::ComputeEvalProofs(batch_digest->value, crs->value);
        if (proofs.size() != out_proof_count) {
            return false;
        }
        for (std::size_t i = 0; i < proofs.size(); ++i) {
            out_proofs[i] = to_boundary_g1(proofs[i]);
        }
        return true;
    } catch (...) {
        return false;
    }
}

bool blst_trx_verify_evaluation_proof(
    const blst_trx_batch_digest* batch_digest,
    const size_t proof_index,
    const blst_p1* proof,
    const blst_trx_crs* crs)
{
    if (batch_digest == nullptr || proof == nullptr || crs == nullptr) {
        return false;
    }

    try {
        return bte::trx::VerifyEvaluationProof(batch_digest->value, proof_index, to_public_g1(*proof), crs->value);
    } catch (...) {
        return false;
    }
}

bool blst_trx_create_partial_decryption(
    const blst_trx_params* params,
    const blst_trx_validator_secret_share* secret_share,
    const blst_trx_batch_digest* batch_digest,
    const blst_trx_public_key_material* public_key_material,
    blst_trx_partial_decryption_share* out_share)
{
    if (params == nullptr || secret_share == nullptr || batch_digest == nullptr ||
        public_key_material == nullptr || out_share == nullptr) {
        return false;
    }

    try {
        const bte::trx::PartialDecryptionShare share = bte::trx::CreatePartialDecryption(
            to_cpp_params(*params),
            to_cpp_secret_share(*secret_share),
            batch_digest->value,
            public_key_material->value);
        out_share->index = share.index;
        out_share->value = to_boundary_g1(share.value);
        return true;
    } catch (...) {
        return false;
    }
}

bool blst_trx_verify_partial_decryption(
    const blst_trx_params* params,
    const blst_trx_partial_decryption_share* share,
    const blst_trx_validator_public_share* public_share,
    const blst_trx_batch_digest* batch_digest,
    const blst_trx_public_key_material* public_key_material)
{
    if (params == nullptr || share == nullptr || public_share == nullptr ||
        batch_digest == nullptr || public_key_material == nullptr) {
        return false;
    }

    try {
        return bte::trx::VerifyPartialDecryption(
            to_cpp_params(*params),
            to_cpp_partial_share(*share),
            to_cpp_public_share(*public_share),
            batch_digest->value,
            public_key_material->value);
    } catch (...) {
        return false;
    }
}

bool blst_trx_combine_shares(
    const blst_trx_partial_decryption_share* shares,
    const size_t share_count,
    const blst_trx_batch_digest* batch_digest,
    const blst_trx_public_key_material* public_key_material,
    blst_p1* out_sigma)
{
    if (shares == nullptr || share_count == 0 || batch_digest == nullptr ||
        public_key_material == nullptr || out_sigma == nullptr) {
        return false;
    }

    try {
        std::vector<bte::trx::PartialDecryptionShare> cpp_shares;
        cpp_shares.reserve(share_count);
        for (std::size_t i = 0; i < share_count; ++i) {
            cpp_shares.push_back(to_cpp_partial_share(shares[i]));
        }
        const blst::G1Point sigma = bte::trx::CombineShares(
            std::span<const bte::trx::PartialDecryptionShare>(cpp_shares.data(), cpp_shares.size()),
            batch_digest->value,
            public_key_material->value);
        *out_sigma = to_boundary_g1(sigma);
        return true;
    } catch (...) {
        return false;
    }
}

bool blst_trx_decapsulate_one(
    const blst_trx_params* params,
    const blst_trx_batch_digest* batch_digest,
    const size_t ciphertext_index,
    const blst_p1* proof,
    const blst_p1* sigma,
    const blst_trx_crs* crs,
    const blst_trx_public_key_material* public_key_material,
    byte* out_session_key,
    const size_t out_session_key_len)
{
    if (params == nullptr || batch_digest == nullptr || proof == nullptr || sigma == nullptr ||
        crs == nullptr || public_key_material == nullptr || out_session_key == nullptr) {
        return false;
    }

    try {
        const auto recovered = bte::trx::DecapsulateOne(
            to_cpp_params(*params),
            batch_digest->value,
            ciphertext_index,
            to_public_g1(*proof),
            to_public_g1(*sigma),
            crs->value,
            public_key_material->value);
        if (!recovered.has_value() || recovered->size() != out_session_key_len) {
            return false;
        }
        std::memcpy(out_session_key, recovered->data(), recovered->size());
        return true;
    } catch (...) {
        return false;
    }
}

void blst_trx_batch_digest_free(blst_trx_batch_digest* batch_digest)
{
    delete batch_digest;
}

}
