#include "BTE/beat/core.h"

#include <array>
#include <cstring>
#include <stdexcept>

#include "BTE/building_blocks/blst_cpp.h"
#include "BTE/building_blocks/hkdf_sha256.h"
#include "BTE/building_blocks/lagrange.h"
#include "blst/fft/middle_product.h"

namespace bte::beat
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

        blst::Fr random_nonzero_fr(const RandomBytesSource& rng, const std::string& dst)
        {
            std::uint32_t counter = 0;
            for (;;) {
                Bytes seed = random_bytes(rng, 48);
                seed.push_back(static_cast<std::uint8_t>(counter >> 24));
                seed.push_back(static_cast<std::uint8_t>(counter >> 16));
                seed.push_back(static_cast<std::uint8_t>(counter >> 8));
                seed.push_back(static_cast<std::uint8_t>(counter));
                ++counter;

                const std::array<std::uint8_t, 32> scalar_bytes =
                    bte::building_blocks::hash_to_scalar_le(seed, dst);
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

        blst::Fr fr_from_u64(const std::uint64_t value)
        {
            return blst::Fr::from_uint64(value);
        }

        blst::Fr fr_add(const blst::Fr& lhs, const blst::Fr& rhs)
        {
            blst::Fr out{};
            blst::Fr::add(out, lhs, rhs);
            return out;
        }

        blst::Fr fr_mul(const blst::Fr& lhs, const blst::Fr& rhs)
        {
            blst::Fr out{};
            blst::Fr::mul(out, lhs, rhs);
            return out;
        }

        blst::G1Point add_g1(const blst::G1Point& lhs, const blst::G1Point& rhs)
        {
            blst::G1Point out{};
            blst::G1Point::dadd(out, lhs, rhs);
            return out;
        }

        blst::GT add_gt(const blst::GT& lhs, const blst::GT& rhs)
        {
            blst::GT out{};
            blst::GT::add(out, lhs, rhs);
            return out;
        }

        blst::GT sub_gt(const blst::GT& lhs, const blst::GT& rhs)
        {
            blst::GT neg_rhs{};
            blst::GT::neg(neg_rhs, rhs);
            return add_gt(lhs, neg_rhs);
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

        Bytes gt_to_bytes(const blst::GT& value)
        {
            std::array<std::uint8_t, 12 * 48> bytes{};
            blst::GT::serialize(bytes.data(), value);
            return Bytes(bytes.begin(), bytes.end());
        }

        Bytes kdf_from_gt(const BeatParameters& params, const blst::GT& value)
        {
            const Bytes ikm = gt_to_bytes(value);
            const std::span<const std::uint8_t> salt(
                reinterpret_cast<const std::uint8_t*>(params.kdf_salt.data()),
                params.kdf_salt.size());
            const std::span<const std::uint8_t> info(
                reinterpret_cast<const std::uint8_t*>(params.kdf_info.data()),
                params.kdf_info.size());
            return bte::building_blocks::hkdf_sha256(salt, ikm, info, params.key_size);
        }

        void validate_setup_params(const BeatParameters& params)
        {
            if (params.max_batch_size == 0) {
                throw std::invalid_argument("BEAT max_batch_size must be nonzero");
            }
            if (params.num_servers == 0) {
                throw std::invalid_argument("BEAT num_servers must be nonzero");
            }
            if (params.threshold == 0 || params.threshold > params.num_servers) {
                throw std::invalid_argument("BEAT threshold must be in [1, num_servers]");
            }
            if (params.key_size == 0) {
                throw std::invalid_argument("BEAT key_size must be nonzero");
            }
        }

        void validate_batch_size(
            const PublicParameters& public_parameters,
            const std::size_t batch_size)
        {
            if (batch_size == 0) {
                throw std::invalid_argument("BEAT batch must be nonempty");
            }
            if (batch_size > public_parameters.max_batch_size) {
                throw std::invalid_argument("BEAT batch exceeds max_batch_size");
            }
        }

        blst::Fr evaluate_polynomial_at(
            const std::span<const blst::Fr> coefficients,
            const blst::Fr& x)
        {
            blst::Fr out = fr_zero();
            for (std::size_t i = coefficients.size(); i-- > 0;) {
                out = fr_add(fr_mul(out, x), coefficients[i]);
            }
            return out;
        }

        blst::G2Point public_g2_power(
            const PublicParameters& public_parameters,
            const std::size_t exponent)
        {
            if (exponent == 0 || exponent >= public_parameters.g2_powers.size()) {
                throw std::invalid_argument("Requested BEAT public power is out of range");
            }
            return public_parameters.g2_powers[exponent];
        }

        blst::G2Point public_g2_centered_offset(
            const PublicParameters& public_parameters,
            const std::ptrdiff_t offset)
        {
            const std::ptrdiff_t min_offset =
                -static_cast<std::ptrdiff_t>(public_parameters.max_batch_size);
            const std::ptrdiff_t max_offset =
                static_cast<std::ptrdiff_t>(public_parameters.max_batch_size) - 1;
            if (offset < min_offset || offset > max_offset) {
                throw std::invalid_argument("Requested BEAT centered CRS offset is out of range");
            }

            const std::size_t exponent = static_cast<std::size_t>(
                static_cast<std::ptrdiff_t>(public_parameters.max_batch_size + 1) + offset);
            return public_g2_power(public_parameters, exponent);
        }

        std::vector<blst::G2Point> make_beta_middle_product_kernel(
            const PublicParameters& public_parameters,
            const std::size_t batch_size)
        {
            std::vector<blst::G2Point> kernel(2 * batch_size, blst::G2Point::infinity());
            for (std::size_t j = 1; j < kernel.size(); ++j) {
                const std::ptrdiff_t offset =
                    static_cast<std::ptrdiff_t>(batch_size) - static_cast<std::ptrdiff_t>(j);
                kernel[j] = public_g2_centered_offset(public_parameters, offset);
            }
            return kernel;
        }

        std::vector<blst::G2Point> make_alpha_middle_product_kernel(
            const PublicParameters& public_parameters,
            const std::size_t batch_size)
        {
            std::vector<blst::G2Point> kernel(2 * batch_size, blst::G2Point::infinity());
            for (std::size_t i = 0; i < batch_size; ++i) {
                kernel[batch_size + i] = public_g2_centered_offset(
                    public_parameters,
                    -static_cast<std::ptrdiff_t>(i + 1));
            }
            return kernel;
        }

        blst::G1Point msm_sum_g1(
            const std::span<const blst::G1Point> points,
            const std::span<const blst::Fr> scalars)
        {
            if (points.size() != scalars.size()) {
                throw std::invalid_argument("BEAT MSM inputs must have matching lengths");
            }
            if (points.empty()) {
                throw std::invalid_argument("BEAT MSM inputs must be nonempty");
            }
            if (points.size() == 1) {
                return scale_g1(points[0], scalars[0]);
            }

            blst::G1Point out = blst::G1Point::infinity();
            for (std::size_t i = 0; i < points.size(); ++i) {
                out = add_g1(out, scale_g1(points[i], scalars[i]));
            }
            return out;
        }

        const blst::MiddleProductPrecomputedB& ensure_precomputed_beta(
            const PreparedPublicParameters& prepared_public_parameters,
            const std::size_t batch_size)
        {
            if (batch_size >= prepared_public_parameters.beta_precomputed_b_by_batch_size.size()) {
                throw std::invalid_argument("BEAT prepared beta cache does not cover this batch size");
            }

            blst::MiddleProductPrecomputedB& slot =
                prepared_public_parameters.beta_precomputed_b_by_batch_size[batch_size];
            if (slot.n != batch_size) {
                const std::vector<blst::G2Point> b = make_beta_middle_product_kernel(
                    prepared_public_parameters.public_parameters,
                    batch_size);
                slot = blst::precompute_middle_product_b(
                    std::span<const blst::G2Point>(b.data(), b.size()),
                    batch_size);
            }
            return slot;
        }

        const blst::MiddleProductPrecomputedB& ensure_precomputed_alpha(
            const PreparedPublicParameters& prepared_public_parameters,
            const std::size_t batch_size)
        {
            if (batch_size >= prepared_public_parameters.alpha_precomputed_b_by_batch_size.size()) {
                throw std::invalid_argument("BEAT prepared alpha cache does not cover this batch size");
            }

            blst::MiddleProductPrecomputedB& slot =
                prepared_public_parameters.alpha_precomputed_b_by_batch_size[batch_size];
            if (slot.n != batch_size) {
                const std::vector<blst::G2Point> b = make_alpha_middle_product_kernel(
                    prepared_public_parameters.public_parameters,
                    batch_size);
                slot = blst::precompute_middle_product_b(
                    std::span<const blst::G2Point>(b.data(), b.size()),
                    batch_size);
            }
            return slot;
        }
    }

    SetupResult Setup(const BeatParameters& params, const RandomBytesSource& rng)
    {
        validate_setup_params(params);

        const blst::Fr tau = random_nonzero_fr(rng, "BTE_BEAT_SETUP_TAU_V1");
        const blst::GT gt_generator = pairing(blst::G1Point::generator(), blst::G2Point::generator());

        std::vector<blst::Fr> tau_powers((2 * params.max_batch_size) + 1, blst::Fr::one());
        for (std::size_t exponent = 1; exponent < tau_powers.size(); ++exponent) {
            tau_powers[exponent] = fr_mul(tau_powers[exponent - 1], tau);
        }

        PublicParameters public_parameters{};
        public_parameters.max_batch_size = params.max_batch_size;
        public_parameters.num_servers = params.num_servers;
        public_parameters.threshold = params.threshold;
        public_parameters.gt_tau_to_bmax_plus_one =
            scale_gt(gt_generator, tau_powers[params.max_batch_size + 1]);
        public_parameters.g2_powers.assign((2 * params.max_batch_size) + 1, blst::G2Point::infinity());
        for (std::size_t exponent = 1; exponent < public_parameters.g2_powers.size(); ++exponent) {
            if (exponent == params.max_batch_size + 1) {
                continue;
            }
            public_parameters.g2_powers[exponent] =
                scale_g2(blst::G2Point::generator(), tau_powers[exponent]);
        }

        std::vector<ServerSecretShare> secret_shares(params.num_servers);
        for (std::size_t server = 0; server < params.num_servers; ++server) {
            secret_shares[server].index = static_cast<std::uint32_t>(server + 1);
            secret_shares[server].tau_power_shares.assign(params.max_batch_size, fr_zero());
        }

        for (std::size_t exponent = 1; exponent <= params.max_batch_size; ++exponent) {
            std::vector<blst::Fr> polynomial(params.threshold, fr_zero());
            polynomial[0] = tau_powers[exponent];
            for (std::size_t degree = 1; degree < params.threshold; ++degree) {
                polynomial[degree] = random_nonzero_fr(rng, "BTE_BEAT_SETUP_SHARE_V1");
            }

            for (std::size_t server = 0; server < params.num_servers; ++server) {
                const blst::Fr x = fr_from_u64(server + 1);
                secret_shares[server].tau_power_shares[exponent - 1] =
                    evaluate_polynomial_at(
                        std::span<const blst::Fr>(polynomial.data(), polynomial.size()),
                        x);
            }
        }

        return SetupResult{ std::move(public_parameters), std::move(secret_shares) };
    }

    PreparedPublicParameters PreparePublicParameters(const PublicParameters& public_parameters)
    {
        PreparedPublicParameters prepared{};
        prepared.public_parameters = public_parameters;
        prepared.beta_precomputed_b_by_batch_size.resize(public_parameters.max_batch_size + 1);
        prepared.alpha_precomputed_b_by_batch_size.resize(public_parameters.max_batch_size + 1);
        return prepared;
    }

    EncryptionResult Encrypt(
        const BeatParameters& params,
        const PublicParameters& public_parameters,
        const RandomBytesSource& rng)
    {
        validate_setup_params(params);
        if (public_parameters.max_batch_size != params.max_batch_size) {
            throw std::invalid_argument("BEAT parameters and public parameters disagree on max_batch_size");
        }

        const blst::Fr r = random_nonzero_fr(rng, "BTE_BEAT_ENCRYPT_R_V1");

        EncryptionResult out{};
        out.ciphertext.c1 = scale_g1(blst::G1Point::generator(), r);
        out.key = kdf_from_gt(params, scale_gt(public_parameters.gt_tau_to_bmax_plus_one, r));
        return out;
    }

    ThresholdRequest BuildThresholdRequest(
        const PublicParameters& public_parameters,
        const std::span<const Ciphertext> ciphertexts)
    {
        validate_batch_size(public_parameters, ciphertexts.size());

        ThresholdRequest request{};
        request.ordered_c1.reserve(ciphertexts.size());
        for (const Ciphertext& ciphertext : ciphertexts) {
            request.ordered_c1.push_back(ciphertext.c1);
        }
        return request;
    }

    PartialDecryptionShare CreatePartialDecryptionShare(
        const ServerSecretShare& secret_share,
        const ThresholdRequest& request)
    {
        if (secret_share.index == 0) {
            throw std::invalid_argument("BEAT server share index must be nonzero");
        }
        if (request.ordered_c1.empty()) {
            throw std::invalid_argument("BEAT threshold request must be nonempty");
        }
        if (secret_share.tau_power_shares.size() < request.ordered_c1.size()) {
            throw std::invalid_argument("BEAT server share does not cover the requested batch size");
        }

        PartialDecryptionShare share{};
        share.index = secret_share.index;
        share.value = msm_sum_g1(
            std::span<const blst::G1Point>(request.ordered_c1.data(), request.ordered_c1.size()),
            std::span<const blst::Fr>(secret_share.tau_power_shares.data(), request.ordered_c1.size()));
        return share;
    }

    blst::G1Point CombinePartialDecryptionShares(
        const std::span<const PartialDecryptionShare> shares,
        const PublicParameters& public_parameters)
    {
        if (shares.empty()) {
            throw std::invalid_argument("Need at least one BEAT share to combine");
        }
        if (shares.size() < public_parameters.threshold) {
            throw std::invalid_argument("Need at least threshold-many BEAT shares to combine");
        }

        std::vector<std::uint32_t> indices;
        indices.reserve(shares.size());
        std::vector<blst::G1Point> share_points;
        share_points.reserve(shares.size());
        for (const PartialDecryptionShare& share : shares) {
            indices.push_back(share.index);
            share_points.push_back(share.value);
        }

        const std::vector<blst::Fr> lambdas = bte::building_blocks::lagrange_at_zero(indices);
        return msm_sum_g1(
            std::span<const blst::G1Point>(share_points.data(), share_points.size()),
            std::span<const blst::Fr>(lambdas.data(), lambdas.size()));
    }

    BatchDecryptionPrecomputation PrecomputeBatchDecryption(
        const PublicParameters& public_parameters,
        const std::span<const Ciphertext> ciphertexts)
    {
        validate_batch_size(public_parameters, ciphertexts.size());

        BatchDecryptionPrecomputation out{};
        out.ciphertexts.assign(ciphertexts.begin(), ciphertexts.end());

        std::vector<blst::G1Point> c1_points;
        c1_points.reserve(ciphertexts.size());
        for (const Ciphertext& ciphertext : ciphertexts) {
            c1_points.push_back(ciphertext.c1);
        }

        std::vector<blst::GT> beta(ciphertexts.size(), blst::GT::identity());
        const std::vector<blst::G2Point> kernel =
            make_beta_middle_product_kernel(public_parameters, ciphertexts.size());
        const auto precomputed_b = blst::precompute_middle_product_b(
            std::span<const blst::G2Point>(kernel.data(), kernel.size()),
            ciphertexts.size());
        blst::MiddleProductWorkspace workspace;
        blst::middle_product(
            std::span<const blst::G1Point>(c1_points.data(), c1_points.size()),
            precomputed_b,
            std::span<blst::GT>(beta.data(), beta.size()),
            workspace);
        out.beta = std::move(beta);
        return out;
    }

    BatchDecryptionPrecomputation PrecomputeBatchDecryption(
        const PreparedPublicParameters& prepared_public_parameters,
        const std::span<const Ciphertext> ciphertexts)
    {
        const PublicParameters& public_parameters = prepared_public_parameters.public_parameters;
        validate_batch_size(public_parameters, ciphertexts.size());

        BatchDecryptionPrecomputation out{};
        out.ciphertexts.assign(ciphertexts.begin(), ciphertexts.end());

        std::vector<blst::G1Point> c1_points;
        c1_points.reserve(ciphertexts.size());
        for (const Ciphertext& ciphertext : ciphertexts) {
            c1_points.push_back(ciphertext.c1);
        }

        std::vector<blst::GT> beta(ciphertexts.size(), blst::GT::identity());
        blst::MiddleProductWorkspace workspace;
        blst::middle_product(
            std::span<const blst::G1Point>(c1_points.data(), c1_points.size()),
            ensure_precomputed_beta(prepared_public_parameters, ciphertexts.size()),
            std::span<blst::GT>(beta.data(), beta.size()),
            workspace);
        out.beta = std::move(beta);
        return out;
    }

    std::vector<Bytes> FinalizeBatchDecryption(
        const BeatParameters& params,
        const PublicParameters& public_parameters,
        const BatchDecryptionPrecomputation& precomputation,
        const blst::G1Point& aggregate_predecryption_key)
    {
        validate_setup_params(params);
        validate_batch_size(public_parameters, precomputation.ciphertexts.size());
        if (precomputation.beta.size() != precomputation.ciphertexts.size()) {
            throw std::invalid_argument("BEAT precomputation beta vector has the wrong size");
        }

        const std::size_t batch_size = precomputation.ciphertexts.size();
        std::vector<blst::GT> cnst(batch_size, blst::GT::identity());
        PreparedPublicParameters prepared{};
        prepared.public_parameters = public_parameters;
        prepared.alpha_precomputed_b_by_batch_size.resize(public_parameters.max_batch_size + 1);
        blst::update(
            aggregate_predecryption_key,
            ensure_precomputed_alpha(prepared, batch_size),
            std::span<blst::GT>(cnst.data(), cnst.size()));

        std::vector<Bytes> keys;
        keys.reserve(batch_size);
        for (std::size_t i = 0; i < batch_size; ++i) {
            keys.push_back(kdf_from_gt(params, sub_gt(cnst[i], precomputation.beta[i])));
        }
        return keys;
    }

    std::vector<Bytes> FinalizeBatchDecryption(
        const BeatParameters& params,
        const PreparedPublicParameters& prepared_public_parameters,
        const BatchDecryptionPrecomputation& precomputation,
        const blst::G1Point& aggregate_predecryption_key)
    {
        validate_setup_params(params);
        validate_batch_size(
            prepared_public_parameters.public_parameters,
            precomputation.ciphertexts.size());
        if (precomputation.beta.size() != precomputation.ciphertexts.size()) {
            throw std::invalid_argument("BEAT precomputation beta vector has the wrong size");
        }

        const std::size_t batch_size = precomputation.ciphertexts.size();
        std::vector<blst::GT> cnst(batch_size, blst::GT::identity());
        blst::update(
            aggregate_predecryption_key,
            ensure_precomputed_alpha(prepared_public_parameters, batch_size),
            std::span<blst::GT>(cnst.data(), cnst.size()));

        std::vector<Bytes> keys;
        keys.reserve(batch_size);
        for (std::size_t i = 0; i < batch_size; ++i) {
            keys.push_back(kdf_from_gt(params, sub_gt(cnst[i], precomputation.beta[i])));
        }
        return keys;
    }

    std::vector<Bytes> RecoverBatchKeys(
        const BeatParameters& params,
        const PublicParameters& public_parameters,
        const std::span<const Ciphertext> ciphertexts,
        const blst::G1Point& aggregate_predecryption_key)
    {
        const BatchDecryptionPrecomputation precomputation =
            PrecomputeBatchDecryption(public_parameters, ciphertexts);
        return FinalizeBatchDecryption(
            params,
            public_parameters,
            precomputation,
            aggregate_predecryption_key);
    }

    std::vector<Bytes> RecoverBatchKeys(
        const BeatParameters& params,
        const PreparedPublicParameters& prepared_public_parameters,
        const std::span<const Ciphertext> ciphertexts,
        const blst::G1Point& aggregate_predecryption_key)
    {
        const BatchDecryptionPrecomputation precomputation =
            PrecomputeBatchDecryption(prepared_public_parameters, ciphertexts);
        return FinalizeBatchDecryption(
            params,
            prepared_public_parameters,
            precomputation,
            aggregate_predecryption_key);
    }
}
