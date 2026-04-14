#include "BTE/trx/core.h"

#include <algorithm>
#include <array>
#include <cstring>
#include <stdexcept>
#include <string>
#include <unordered_set>

#include "BTE/building_blocks/bls_ots.h"
#include "BTE/building_blocks/blst_cpp.h"
#include "BTE/building_blocks/hkdf_sha256.h"
#include "BTE/building_blocks/kzg.h"
#include "BTE/building_blocks/lagrange.h"
#include "BTE/building_blocks/polynomial.h"
#include "BTE/building_blocks/serialization.h"

namespace bte::trx
{
    namespace
    {
        blst::Fr fr_zero()
        {
            blst::Fr zero{};
            std::memset(&zero, 0, sizeof(zero));
            return zero;
        }

        bool fr_is_zero(const blst::Fr& value)
        {
            return blst::Fr::equal(value, fr_zero());
        }

        blst::Fr bytes_to_fr(const std::array<std::uint8_t, 32>& bytes)
        {
            blst::Fr out{};
            blst::Fr::from_scalar(out, bytes.data());
            return out;
        }

        Bytes random_bytes(const RandomBytesSource& rng, const std::size_t count)
        {
            if (!rng) {
                throw std::invalid_argument("Random source is required");
            }

            Bytes out(count);
            if (!out.empty()) {
                rng(std::span<std::uint8_t>(out.data(), out.size()));
            }
            return out;
        }

        blst::Fr random_nonzero_fr(const RandomBytesSource& rng)
        {
            if (!rng) {
                throw std::invalid_argument("Random source is required");
            }

            std::uint32_t counter = 0;
            for (;;) {
                Bytes seed = random_bytes(rng, 48);
                bte::building_blocks::append_u32_be(seed, counter++);
                const std::array<std::uint8_t, 32> scalar_bytes =
                    bte::building_blocks::hash_to_scalar_le(seed, "BTE_TRX_RNG_FR_V1");
                const blst::Fr out = bytes_to_fr(scalar_bytes);
                if (!fr_is_zero(out)) {
                    return out;
                }
            }
        }

        std::array<std::uint8_t, 32> fr_to_scalar_bytes(const blst::Fr& scalar)
        {
            std::array<std::uint8_t, 32> bytes{};
            blst::Fr::to_scaler(bytes.data(), scalar);
            return bytes;
        }

        blst::G1Point add_g1(const blst::G1Point& lhs, const blst::G1Point& rhs)
        {
            blst::G1Point out{};
            blst::G1Point::dadd(out, lhs, rhs);
            return out;
        }

        blst::G2Point add_g2(const blst::G2Point& lhs, const blst::G2Point& rhs)
        {
            blst::G2Point out{};
            blst::G2Point::dadd(out, lhs, rhs);
            return out;
        }

        blst::G1Point sub_g1(const blst::G1Point& lhs, const blst::G1Point& rhs)
        {
            blst::G1Point neg_rhs{};
            blst::G1Point::neg(neg_rhs, rhs);
            return add_g1(lhs, neg_rhs);
        }

        blst::G2Point sub_g2(const blst::G2Point& lhs, const blst::G2Point& rhs)
        {
            blst::G2Point neg_rhs{};
            blst::G2Point::neg(neg_rhs, rhs);
            return add_g2(lhs, neg_rhs);
        }

        blst::G1Point scale_g1(const blst::G1Point& point, const blst::Fr& scalar)
        {
            if (fr_is_zero(scalar)) {
                return blst::G1Point::infinity();
            }

            const std::array<std::uint8_t, 32> scalar_bytes = fr_to_scalar_bytes(scalar);
            blst::G1Point out{};
            blst::G1Point::mult(out, point, scalar_bytes.data(), 256);
            return out;
        }

        blst::G2Point scale_g2(const blst::G2Point& point, const blst::Fr& scalar)
        {
            if (fr_is_zero(scalar)) {
                return blst::G2Point::infinity();
            }

            const std::array<std::uint8_t, 32> scalar_bytes = fr_to_scalar_bytes(scalar);
            blst::G2Point out{};
            blst::G2Point::mult(out, point, scalar_bytes.data(), 256);
            return out;
        }

        blst::GT scale_gt(const blst::GT& value, const blst::Fr& scalar)
        {
            if (fr_is_zero(scalar)) {
                return blst::GT::identity();
            }

            const std::array<std::uint8_t, 32> scalar_bytes = fr_to_scalar_bytes(scalar);
            blst::GT out{};
            blst::GT::mult(out, value, scalar_bytes.data(), 256);
            return out;
        }

        blst::Fr fr_mul(const blst::Fr& lhs, const blst::Fr& rhs)
        {
            blst::Fr out{};
            blst::Fr::mul(out, lhs, rhs);
            return out;
        }

        blst::Fr fr_neg(const blst::Fr& value)
        {
            blst::Fr out{};
            blst::Fr::neg(out, value);
            return out;
        }

        blst::GT pairing(const blst::G1Point& p, const blst::G2Point& q)
        {
            blst::G1AffinePoint p_aff{};
            blst::G2AffinePoint q_aff{};
            blst::G1Point::to_affine(p_aff, p);
            blst::G2Point::to_affine(q_aff, q);

            blst::GT miller{};
            blst::GT::miller_loop(miller, &q_aff, &p_aff);
            blst::GT out{};
            blst::GT::final_mult(out, miller);
            return out;
        }

        std::array<std::uint8_t, 48> compress_g1(const blst::G1Point& point)
        {
            std::array<std::uint8_t, 48> out{};
            blst::G1Point::compress(out.data(), point);
            return out;
        }

        std::array<std::uint8_t, 96> compress_g2(const blst::G2Point& point)
        {
            std::array<std::uint8_t, 96> out{};
            blst::G2Point::compress(out.data(), point);
            return out;
        }

        blst::G1Point deserialize_g1(const std::array<std::uint8_t, 48>& bytes)
        {
            blst::G1AffinePoint affine{};
            if (blst::G1AffinePoint::uncompress(affine, bytes.data()) != blst::BLST_SUCCESS) {
                throw std::invalid_argument("Invalid compressed G1 point");
            }
            if (blst::G1AffinePoint::is_inf(affine) != 0 || blst::G1AffinePoint::in_subgroup(affine) == 0) {
                throw std::invalid_argument("Compressed G1 point is not a valid subgroup element");
            }
            blst::G1Point out{};
            blst::G1Point::from_affine(out, affine);
            return out;
        }

        blst::G2Point deserialize_g2(const std::array<std::uint8_t, 96>& bytes)
        {
            blst::G2AffinePoint affine{};
            if (blst::G2AffinePoint::uncompress(affine, bytes.data()) != blst::BLST_SUCCESS) {
                throw std::invalid_argument("Invalid compressed G2 point");
            }
            if (blst::G2AffinePoint::is_inf(affine) != 0 || blst::G2AffinePoint::in_subgroup(affine) == 0) {
                throw std::invalid_argument("Compressed G2 point is not a valid subgroup element");
            }
            blst::G2Point out{};
            blst::G2Point::from_affine(out, affine);
            return out;
        }

        Bytes serialize_ciphertext_header(
            const TrxParameters& params,
            const std::span<const std::uint8_t> associated_data,
            const std::array<std::uint8_t, 96>& ct1,
            const std::array<std::uint8_t, 96>& ct2,
            const Bytes& ct3,
            const std::array<std::uint8_t, 48>& vk_sig)
        {
            Bytes out;
            out.reserve(1 + 4 + associated_data.size() + ct1.size() + ct2.size() + 4 + ct3.size() + vk_sig.size());
            bte::building_blocks::append_u8(out, params.serialization_version);
            bte::building_blocks::append_u32_be(out, static_cast<std::uint32_t>(associated_data.size()));
            bte::building_blocks::append_bytes(out, associated_data);
            bte::building_blocks::append_array(out, ct1);
            bte::building_blocks::append_array(out, ct2);
            bte::building_blocks::append_u32_be(out, static_cast<std::uint32_t>(ct3.size()));
            bte::building_blocks::append_bytes(out, std::span<const std::uint8_t>(ct3.data(), ct3.size()));
            bte::building_blocks::append_array(out, vk_sig);
            return out;
        }

        Bytes serialize_tag_input(
            const std::array<std::uint8_t, 48>& vk_sig,
            const std::span<const std::uint8_t> associated_data)
        {
            Bytes out;
            out.reserve(vk_sig.size() + 4 + associated_data.size());
            bte::building_blocks::append_array(out, vk_sig);
            bte::building_blocks::append_u32_be(out, static_cast<std::uint32_t>(associated_data.size()));
            bte::building_blocks::append_bytes(out, associated_data);
            return out;
        }

        Bytes g2_to_bytes(const blst::G2Point& point)
        {
            const std::array<std::uint8_t, 96> compressed = compress_g2(point);
            return Bytes(compressed.begin(), compressed.end());
        }

        Bytes gt_to_bytes(const blst::GT& value)
        {
            std::array<std::uint8_t, 12 * 48> bytes{};
            blst::GT::serialize(bytes.data(), value);
            return Bytes(bytes.begin(), bytes.end());
        }

        blst::Fr hash_to_fr(const Bytes& message, const std::string& dst)
        {
            return bytes_to_fr(bte::building_blocks::hash_to_scalar_le(message, dst));
        }

        blst::G1Point hash_to_g1(const Bytes& message, const std::string& dst)
        {
            return bte::building_blocks::hash_to_g1(message, dst);
        }

        template<typename KeyLike>
        const blst::G1Point& cached_h1_for_pk(
            const KeyLike& key_material,
            const TrxParameters& params)
        {
            if (!key_material.has_cached_h1_pk || key_material.cached_h1_dst != params.h1_dst) {
                key_material.cached_h1_pk = hash_to_g1(g2_to_bytes(key_material.pk), params.h1_dst);
                key_material.cached_h1_dst = params.h1_dst;
                key_material.has_cached_h1_pk = true;
            }
            return key_material.cached_h1_pk;
        }

        Bytes kdf_from_gt(const TrxParameters& params, const blst::GT& value, const std::size_t output_len)
        {
            const Bytes gt_bytes = gt_to_bytes(value);
            const Bytes salt(params.kdf_salt.begin(), params.kdf_salt.end());
            const Bytes info(params.kdf_info.begin(), params.kdf_info.end());
            return bte::building_blocks::hkdf_sha256(
                std::span<const std::uint8_t>(salt.data(), salt.size()),
                std::span<const std::uint8_t>(gt_bytes.data(), gt_bytes.size()),
                std::span<const std::uint8_t>(info.data(), info.size()),
                output_len);
        }

        Bytes xor_bytes(const Bytes& lhs, const Bytes& rhs)
        {
            if (lhs.size() != rhs.size()) {
                throw std::invalid_argument("XOR operands must have the same length");
            }

            Bytes out(lhs.size());
            for (std::size_t i = 0; i < lhs.size(); ++i) {
                out[i] = lhs[i] ^ rhs[i];
            }
            return out;
        }

        std::vector<ValidatorPublicShare> collect_subset(
            const std::span<const ValidatorPublicShare> shares,
            const std::span<const std::uint32_t> subset_indices)
        {
            if (shares.empty()) {
                throw std::invalid_argument("Need at least one public share");
            }

            if (subset_indices.empty()) {
                return std::vector<ValidatorPublicShare>(shares.begin(), shares.end());
            }

            std::vector<ValidatorPublicShare> subset;
            subset.reserve(subset_indices.size());
            for (const std::uint32_t index : subset_indices) {
                const auto it = std::find_if(
                    shares.begin(),
                    shares.end(),
                    [index](const ValidatorPublicShare& share) { return share.index == index; });
                if (it == shares.end()) {
                    throw std::invalid_argument("Subset index not present in public shares");
                }
                subset.push_back(*it);
            }
            return subset;
        }

        blst::GT derive_ciphertext_secret(
            const VerifiedCiphertext& ciphertext,
            const blst::G1Point& proof,
            const blst::G1Point& sigma)
        {
            const blst::GT left = pairing(proof, ciphertext.ct1);
            const blst::GT right = pairing(sigma, ciphertext.ct2);
            blst::GT out{};
            blst::GT::add(out, left, right);
            return out;
        }
    }

    Crs GenerateCRS(const TrxParameters& params, const RandomBytesSource& rng)
    {
        if (params.max_batch_size == 0 || params.num_contexts == 0) {
            throw std::invalid_argument("TrX parameters must define a nonzero batch size and context count");
        }

        const blst::Fr tau = random_nonzero_fr(rng);

        Crs crs{};
        crs.max_batch_size = params.max_batch_size;
        crs.contexts.resize(params.num_contexts);
        crs.h_tau = scale_g2(blst::G2Point::generator(), tau);

        for (std::size_t context_id = 0; context_id < params.num_contexts; ++context_id) {
            const blst::Fr kappa = random_nonzero_fr(rng);
            blst::Fr tau_power = blst::Fr::one();

            CrsContext context{};
            context.g1_powers.reserve(params.max_batch_size + 1);
            for (std::size_t degree = 0; degree <= params.max_batch_size; ++degree) {
                context.g1_powers.push_back(
                    scale_g1(blst::G1Point::generator(), fr_mul(kappa, tau_power)));
                tau_power = fr_mul(tau_power, tau);
            }
            crs.contexts[context_id] = std::move(context);
        }

        return crs;
    }

    ValidatorPublicShare DeriveValidatorPublicShare(
        const ValidatorSecretShare& secret_share,
        const blst::G2Point& h_tau)
    {
        if (secret_share.index == 0) {
            throw std::invalid_argument("Validator index must be nonzero");
        }

        ValidatorPublicShare public_share{};
        public_share.index = secret_share.index;
        public_share.pk_share = scale_g2(blst::G2Point::generator(), secret_share.sk_share);
        public_share.pk_tau_share = scale_g2(h_tau, secret_share.sk_share);
        return public_share;
    }

    PublicKeyMaterial DerivePublicKeyMaterial(
        const std::span<const ValidatorPublicShare> public_shares,
        const std::span<const std::uint32_t> subset_indices)
    {
        const std::vector<ValidatorPublicShare> subset = collect_subset(public_shares, subset_indices);
        std::vector<std::uint32_t> indices;
        indices.reserve(subset.size());
        for (const ValidatorPublicShare& share : subset) {
            indices.push_back(share.index);
        }

        const std::vector<blst::Fr> lambdas = bte::building_blocks::lagrange_at_zero(indices);

        PublicKeyMaterial material{};
        material.validator_public_shares.assign(public_shares.begin(), public_shares.end());
        material.pk = blst::G2Point::infinity();
        material.pk_tau = blst::G2Point::infinity();
        material.threshold = subset_indices.empty() ? public_shares.size() : subset.size();
        for (std::size_t i = 0; i < subset.size(); ++i) {
            material.pk = add_g2(material.pk, scale_g2(subset[i].pk_share, lambdas[i]));
            material.pk_tau = add_g2(material.pk_tau, scale_g2(subset[i].pk_tau_share, lambdas[i]));
        }
        return material;
    }

    EncryptionKey ToEncryptionKey(const PublicKeyMaterial& public_key_material)
    {
        EncryptionKey key{};
        key.pk = public_key_material.pk;
        key.pk_tau = public_key_material.pk_tau;
        return key;
    }

    EncryptedCiphertext Encrypt(
        const TrxParameters& params,
        const EncryptionKey& encryption_key,
        const std::span<const std::uint8_t> associated_data,
        const RandomBytesSource& rng)
    {
        if (params.session_key_size == 0) {
            throw std::invalid_argument("Session key size must be nonzero");
        }

        EncryptedCiphertext output{};
        output.session_key = random_bytes(rng, params.session_key_size);

        const Bytes one_time_ikm = random_bytes(rng, 32);
        const Bytes placeholder_message{};
        const bte::building_blocks::BlsMinPkOneTimeSignature placeholder =
            bte::building_blocks::sign_bls_min_pk(
                std::span<const std::uint8_t>(one_time_ikm.data(), one_time_ikm.size()),
                std::span<const std::uint8_t>(placeholder_message.data(), placeholder_message.size()),
                params.ots_sign_dst);

        output.associated_data.assign(associated_data.begin(), associated_data.end());
        const blst::Fr tag = hash_to_fr(
            serialize_tag_input(placeholder.verification_key, associated_data),
            params.hf_dst);
        const blst::Fr alpha = random_nonzero_fr(rng);

        const blst::G2Point neg_tag_pk = scale_g2(encryption_key.pk, fr_neg(tag));
        const blst::G2Point ct1_base = add_g2(encryption_key.pk_tau, neg_tag_pk);
        const blst::G2Point ct1_point = scale_g2(ct1_base, alpha);
        const blst::G2Point ct2_point = scale_g2(blst::G2Point::generator(), alpha);

        const blst::G1Point& h1 = cached_h1_for_pk(encryption_key, params);
        const blst::GT base_secret = pairing(h1, encryption_key.pk);
        const blst::GT encapsulation_secret = scale_gt(base_secret, alpha);
        const Bytes mask = kdf_from_gt(params, encapsulation_secret, output.session_key.size());
        const Bytes ct3 = xor_bytes(output.session_key, mask);

        output.ciphertext.ct1 = compress_g2(ct1_point);
        output.ciphertext.ct2 = compress_g2(ct2_point);
        output.ciphertext.ct3 = ct3;
        output.ciphertext.vk_sig = placeholder.verification_key;

        const Bytes header = serialize_ciphertext_header(
            params,
            associated_data,
            output.ciphertext.ct1,
            output.ciphertext.ct2,
            output.ciphertext.ct3,
            output.ciphertext.vk_sig);
        const bte::building_blocks::BlsMinPkOneTimeSignature signed_header =
            bte::building_blocks::sign_bls_min_pk(
                std::span<const std::uint8_t>(one_time_ikm.data(), one_time_ikm.size()),
                std::span<const std::uint8_t>(header.data(), header.size()),
                params.ots_sign_dst);
        output.ciphertext.sig = signed_header.signature;
        return output;
    }

    EncryptedCiphertext Encrypt(
        const TrxParameters& params,
        const EncryptionKey& encryption_key,
        const RandomBytesSource& rng)
    {
        return Encrypt(params, encryption_key, std::span<const std::uint8_t>(), rng);
    }

    std::optional<VerifiedCiphertext> VerifyCiphertext(
        const TrxParameters& params,
        const Ciphertext& ciphertext,
        const std::span<const std::uint8_t> associated_data)
    {
        try {
            const blst::G2Point ct1 = deserialize_g2(ciphertext.ct1);
            const blst::G2Point ct2 = deserialize_g2(ciphertext.ct2);

            const Bytes header = serialize_ciphertext_header(
                params,
                associated_data,
                ciphertext.ct1,
                ciphertext.ct2,
                ciphertext.ct3,
                ciphertext.vk_sig);

            if (!bte::building_blocks::verify_bls_min_pk(
                    std::span<const std::uint8_t>(ciphertext.vk_sig.data(), ciphertext.vk_sig.size()),
                    std::span<const std::uint8_t>(ciphertext.sig.data(), ciphertext.sig.size()),
                    std::span<const std::uint8_t>(header.data(), header.size()),
                    params.ots_sign_dst)) {
                return std::nullopt;
            }

            VerifiedCiphertext verified{};
            verified.ciphertext = ciphertext;
            verified.associated_data.assign(associated_data.begin(), associated_data.end());
            verified.ct1 = ct1;
            verified.ct2 = ct2;
            verified.tag = hash_to_fr(
                serialize_tag_input(ciphertext.vk_sig, associated_data),
                params.hf_dst);
            return verified;
        } catch (...) {
            return std::nullopt;
        }
    }

    std::optional<VerifiedCiphertext> VerifyCiphertext(
        const TrxParameters& params,
        const Ciphertext& ciphertext)
    {
        return VerifyCiphertext(params, ciphertext, std::span<const std::uint8_t>());
    }

    BatchDigest BuildBatchDigest(
        const TrxParameters& params,
        const std::span<const CiphertextEnvelope> ciphertexts,
        const std::size_t context_id,
        const Crs& crs)
    {
        if (context_id >= crs.contexts.size()) {
            throw std::invalid_argument("Context id is out of range");
        }
        if (ciphertexts.empty()) {
            throw std::invalid_argument("Need at least one ciphertext");
        }

        BatchDigest digest{};
        digest.context_id = context_id;

        for (std::size_t i = 0; i < ciphertexts.size(); ++i) {
            const std::optional<VerifiedCiphertext> verified = VerifyCiphertext(
                params,
                ciphertexts[i].ciphertext,
                std::span<const std::uint8_t>(
                    ciphertexts[i].associated_data.data(),
                    ciphertexts[i].associated_data.size()));
            if (!verified.has_value()) {
                continue;
            }
            digest.valid_indices.push_back(i);
            digest.valid_ciphertexts.push_back(*verified);
            digest.tags.push_back(verified->tag);
        }

        if (digest.valid_ciphertexts.empty()) {
            throw std::runtime_error("No valid ciphertexts remain after verification");
        }
        if (digest.valid_ciphertexts.size() > params.max_batch_size) {
            throw std::runtime_error("Batch exceeds configured maximum size");
        }

        std::unordered_set<std::string> seen_tags;
        seen_tags.reserve(digest.tags.size());
        for (const blst::Fr& tag : digest.tags) {
            const std::array<std::uint8_t, 32> tag_bytes = fr_to_scalar_bytes(tag);
            const std::string key(tag_bytes.begin(), tag_bytes.end());
            if (!seen_tags.insert(key).second) {
                throw std::runtime_error("Duplicate tags are rejected in v1");
            }
        }

        digest.polynomial = bte::building_blocks::build_vanishing_polynomial(digest.tags);
        digest.commitment = bte::building_blocks::commit_polynomial(
            std::span<const blst::G1Point>(
                crs.contexts[context_id].g1_powers.data(),
                crs.contexts[context_id].g1_powers.size()),
            digest.polynomial);
        return digest;
    }

    BatchDigest BuildBatchDigest(
        const TrxParameters& params,
        const std::span<const Ciphertext> ciphertexts,
        const std::size_t context_id,
        const Crs& crs)
    {
        std::vector<CiphertextEnvelope> envelopes;
        envelopes.reserve(ciphertexts.size());
        for (const Ciphertext& ciphertext : ciphertexts) {
            envelopes.push_back({ ciphertext, {} });
        }
        return BuildBatchDigest(
            params,
            std::span<const CiphertextEnvelope>(envelopes.data(), envelopes.size()),
            context_id,
            crs);
    }

    std::vector<blst::G1Point> ComputeEvalProofs(
        const BatchDigest& digest,
        const Crs& crs)
    {
        if (digest.context_id >= crs.contexts.size()) {
            throw std::invalid_argument("Context id is out of range");
        }
        if (digest.has_cached_eval_proofs) {
            return digest.cached_eval_proofs;
        }

        const auto& powers = crs.contexts[digest.context_id].g1_powers;
        std::vector<blst::G1Point> proofs = bte::building_blocks::compute_fk_eval_proofs(
            std::span<const blst::G1Point>(powers.data(), powers.size()),
            std::span<const blst::Fr>(digest.tags.data(), digest.tags.size()));
        digest.cached_eval_proofs = proofs;
        digest.has_cached_eval_proofs = true;
        return proofs;
    }

    bool VerifyEvaluationProof(
        const BatchDigest& digest,
        const std::size_t proof_index,
        const blst::G1Point& proof,
        const Crs& crs)
    {
        if (proof_index >= digest.tags.size()) {
            throw std::invalid_argument("Proof index out of range");
        }
        return bte::building_blocks::verify_zero_opening(
            digest.commitment,
            proof,
            digest.tags[proof_index],
            crs.h_tau);
    }

    PartialDecryptionShare CreatePartialDecryption(
        const TrxParameters& params,
        const ValidatorSecretShare& secret_share,
        const BatchDigest& digest,
        const PublicKeyMaterial& public_key_material)
    {
        if (secret_share.index == 0) {
            throw std::invalid_argument("Validator index must be nonzero");
        }

        PartialDecryptionShare share{};
        share.index = secret_share.index;

        const blst::G1Point& h1 = cached_h1_for_pk(public_key_material, params);
        const blst::G1Point base = sub_g1(h1, digest.commitment);
        share.value = scale_g1(base, secret_share.sk_share);
        return share;
    }

    bool VerifyPartialDecryption(
        const TrxParameters& params,
        const PartialDecryptionShare& share,
        const ValidatorPublicShare& public_share,
        const BatchDigest& digest,
        const PublicKeyMaterial& public_key_material)
    {
        if (share.index != public_share.index) {
            return false;
        }

        const blst::G1Point& h1 = cached_h1_for_pk(public_key_material, params);
        const blst::G1Point base = sub_g1(h1, digest.commitment);

        const blst::GT lhs = pairing(base, public_share.pk_share);
        const blst::GT rhs = pairing(share.value, blst::G2Point::generator());
        return blst::GT::equal(lhs, rhs) != 0;
    }

    blst::G1Point CombineShares(
        const std::span<const PartialDecryptionShare> shares,
        const BatchDigest& digest,
        const PublicKeyMaterial& public_key_material)
    {
        (void)digest;
        (void)public_key_material;

        if (shares.empty()) {
            throw std::invalid_argument("Need at least one share to combine");
        }
        if (public_key_material.threshold == 0) {
            throw std::invalid_argument("Public key material must specify a nonzero threshold");
        }
        if (shares.size() < public_key_material.threshold) {
            throw std::invalid_argument("Need at least threshold-many shares to combine");
        }

        std::vector<std::uint32_t> indices;
        indices.reserve(shares.size());
        for (const PartialDecryptionShare& share : shares) {
            indices.push_back(share.index);
        }
        const std::vector<blst::Fr> lambdas = bte::building_blocks::lagrange_at_zero(indices);

        blst::G1Point sigma = blst::G1Point::infinity();
        for (std::size_t i = 0; i < shares.size(); ++i) {
            sigma = add_g1(sigma, scale_g1(shares[i].value, lambdas[i]));
        }
        return sigma;
    }

    std::optional<Bytes> DecapsulateOne(
        const TrxParameters& params,
        const BatchDigest& digest,
        const std::size_t ciphertext_index,
        const blst::G1Point& proof,
        const blst::G1Point& sigma,
        const Crs& crs,
        const PublicKeyMaterial& public_key_material)
    {
        (void)public_key_material;
        if (ciphertext_index >= digest.valid_ciphertexts.size()) {
            throw std::invalid_argument("Ciphertext index is out of range");
        }
        if (!VerifyEvaluationProof(digest, ciphertext_index, proof, crs)) {
            return std::nullopt;
        }

        const VerifiedCiphertext& ciphertext = digest.valid_ciphertexts[ciphertext_index];
        const blst::GT secret = derive_ciphertext_secret(ciphertext, proof, sigma);
        const Bytes mask = kdf_from_gt(params, secret, ciphertext.ciphertext.ct3.size());
        return xor_bytes(ciphertext.ciphertext.ct3, mask);
    }
}
