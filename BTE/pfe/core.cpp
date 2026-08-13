#include "BTE/pfe/core.h"

#include <array>
#include <cstring>
#include <stdexcept>

#include "BTE/building_blocks/blst_cpp.h"
#include "BTE/building_blocks/hkdf_sha256.h"
#include "BTE/building_blocks/lagrange.h"
#include "blst/fft/context.h"
#include "blst/fft/ctFft.h"
#include "blst/fft/primitive_root.h"
#if BLST_DEFAULT_AVX
#include "blst/avx/curve/map_to_G.h"
#endif
#include "blst/curve/multi_scalar.h"

namespace bte::pfe
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

        blst::Fr fr_from_i64(const std::int64_t value)
        {
            if (value >= 0) {
                return fr_from_u64(static_cast<std::uint64_t>(value));
            }

            blst::Fr out = fr_from_u64(static_cast<std::uint64_t>(-value));
            blst::Fr::neg(out, out);
            return out;
        }

        blst::Fr fr_add(const blst::Fr& lhs, const blst::Fr& rhs)
        {
            blst::Fr out{};
            blst::Fr::add(out, lhs, rhs);
            return out;
        }

        blst::Fr fr_sub(const blst::Fr& lhs, const blst::Fr& rhs)
        {
            blst::Fr out{};
            blst::Fr::sub(out, lhs, rhs);
            return out;
        }

        blst::Fr fr_mul(const blst::Fr& lhs, const blst::Fr& rhs)
        {
            blst::Fr out{};
            blst::Fr::mul(out, lhs, rhs);
            return out;
        }

        blst::Fr fr_inverse(const blst::Fr& value)
        {
            if (fr_is_zero(value)) {
                throw std::invalid_argument("Cannot invert zero");
            }

            blst::Fr out{};
            blst::Fr::reciprocal(out, value);
            return out;
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

        blst::GT add_gt(const blst::GT& lhs, const blst::GT& rhs)
        {
            blst::GT out{};
            blst::GT::add(out, lhs, rhs);
            return out;
        }

        blst::G1Point sub_g1(const blst::G1Point& lhs, const blst::G1Point& rhs)
        {
            blst::G1Point neg_rhs{};
            blst::G1Point::neg(neg_rhs, rhs);
            return add_g1(lhs, neg_rhs);
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

        blst::G1Point msm_sum_g1(
            const std::span<const blst::G1Point> points,
            const std::span<const blst::Fr> scalars)
        {
            if (points.size() != scalars.size()) {
                throw std::invalid_argument("PFE MSM inputs must have matching lengths");
            }
            if (points.empty()) {
                throw std::invalid_argument("PFE MSM inputs must be nonempty");
            }
            if (points.size() == 1) {
                return scale_g1(points[0], scalars[0]);
            }
            if (points.size() < 8) {
                blst::G1Point out = blst::G1Point::infinity();
                for (std::size_t i = 0; i < points.size(); ++i) {
                    out = add_g1(out, scale_g1(points[i], scalars[i]));
                }
                return out;
            }

            std::vector<blst::scalar::G1AffinePoint> affine_points(points.size());
            std::vector<std::uint8_t> scalar_blob(points.size() * 32, 0);
            for (std::size_t i = 0; i < points.size(); ++i) {
#if BLST_DEFAULT_AVX
                const blst::scalar::G1Point scalar_point =
                    blst::avx::avx_map_to_G_detail::to_scalar(points[i]);
                blst::scalar::G1Point::to_affine(affine_points[i], scalar_point);
#else
                blst::scalar::G1Point::to_affine(affine_points[i], points[i]);
#endif
                const std::array<std::uint8_t, 32> scalar_bytes = fr_to_scalar_bytes(scalars[i]);
                std::memcpy(scalar_blob.data() + (i * 32), scalar_bytes.data(), scalar_bytes.size());
            }

            blst::scalar::G1AffinePoint* point_params[2] = { affine_points.data(), nullptr };
            const blst::byte* scalar_params[2] = { scalar_blob.data(), nullptr };

            const std::size_t scratch_bytes =
                blst::msm_mult_pippenger_scratch_sizeof_<blst::Pointxyzz<blst::Fp>>(points.size());
            std::vector<blst::Pointxyzz<blst::Fp>> scratch(
                scratch_bytes / sizeof(blst::Pointxyzz<blst::Fp>));

            blst::scalar::G1Point scalar_out{};
            blst::msm_mult_pippenger<
                blst::scalar::Point<blst::scalar::Fp>,
                blst::scalar::AffinePoint<blst::scalar::Fp>,
                blst::scalar::Pointxyzz<blst::scalar::Fp>>(
                scalar_out,
                point_params,
                points.size(),
                scalar_params,
                256,
                scratch.data(),
                blst::pippenger_window_size(points.size()));
#if BLST_DEFAULT_AVX
            return blst::avx::avx_map_to_G_detail::from_scalar(scalar_out);
#else
            return scalar_out;
#endif
        }

        blst::GT pairing(const blst::G1Point& p, const blst::G2Point& q)
        {
            blst::G1AffinePoint p_aff{};
            blst::G2AffinePoint q_aff{};
            blst::G1Point::to_affine(p_aff, p);
            blst::G2Point::to_affine(q_aff, q);

            blst::MillerResult miller{};
            blst::GT::miller_loop(miller, &q_aff, &p_aff);
            blst::GT out{};
            blst::GT::final_exp(out, miller);
            return out;
        }

        Bytes gt_to_bytes(const blst::GT& value)
        {
            std::array<std::uint8_t, 12 * 48> bytes{};
            blst::GT::serialize(bytes.data(), value);
            return Bytes(bytes.begin(), bytes.end());
        }

        Bytes kdf_from_gt(const PfeParameters& params, const blst::GT& value)
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

        blst::Fr evaluate_polynomial_at(
            const std::span<const blst::Fr> coefficients,
            const blst::Fr& point)
        {
            blst::Fr out = fr_zero();
            for (auto it = coefficients.rbegin(); it != coefficients.rend(); ++it) {
                out = fr_add(fr_mul(out, point), *it);
            }
            return out;
        }

        void validate_params(const PfeParameters& params)
        {
            if (params.batch_size == 0) {
                throw std::invalid_argument("PFE batch_size must be nonzero");
            }
            if (params.num_servers == 0) {
                throw std::invalid_argument("PFE num_servers must be nonzero");
            }
            if (params.threshold == 0 || params.threshold > params.num_servers) {
                throw std::invalid_argument("PFE threshold must be in [1, num_servers]");
            }
            if (params.key_size == 0) {
                throw std::invalid_argument("PFE key_size must be nonzero");
            }
        }

        void validate_batch_size(
            const PublicParameters& public_parameters,
            const std::size_t batch_size)
        {
            if (batch_size != public_parameters.batch_size) {
                throw std::invalid_argument("PFE batches must match the fixed setup batch size");
            }
        }

        blst::Fr partial_fraction_value(const blst::Fr& x, const std::int64_t i)
        {
            return fr_inverse(fr_add(x, fr_from_i64(i)));
        }

        bool sample_is_valid_setup_secret(const blst::Fr& x, const std::size_t batch_size)
        {
            if (fr_is_zero(x)) {
                return false;
            }
            if (fr_is_zero(fr_sub(x, blst::Fr::one()))) {
                return false;
            }
            for (std::size_t i = 1; i <= batch_size; ++i) {
                if (fr_is_zero(fr_add(x, fr_from_u64(i)))) {
                    return false;
                }
            }
            return true;
        }

        blst::Fr sample_valid_setup_secret(
            const PfeParameters& params,
            const RandomBytesSource& rng)
        {
            for (;;) {
                const blst::Fr candidate = random_nonzero_fr(rng, "BTE_PFE_SETUP_X_V1");
                if (sample_is_valid_setup_secret(candidate, params.batch_size)) {
                    return candidate;
                }
            }
        }

        std::size_t next_power_of_two(std::size_t value)
        {
            std::size_t out = 1;
            while (out < value) {
                out <<= 1;
            }
            return out;
        }

        std::vector<blst::Fr> make_fft_roots(const std::size_t n)
        {
            const blst::Fr primitive = blst::primRootOfUnity<blst::Fr>(n);
            std::vector<blst::Fr> roots(n, fr_zero());
            blst::Fr current = blst::Fr::one();
            for (std::size_t i = 0; i < n; ++i) {
                roots[i] = current;
                blst::Fr::mul(current, current, primitive);
            }
            return roots;
        }

        std::vector<blst::Fr> make_ifft_roots(const std::span<const blst::Fr> roots)
        {
            std::vector<blst::Fr> inv_roots(roots.size(), fr_zero());
            for (std::size_t i = 0; i < roots.size(); ++i) {
                blst::Fr::reciprocal(inv_roots[i], roots[i]);
            }
            return inv_roots;
        }

        std::vector<blst::Fr> build_toeplitz_kernel_embedding(
            const std::size_t batch_size,
            const std::size_t fft_size)
        {
            std::vector<blst::Fr> kernel(fft_size, fr_zero());
            for (std::size_t i = 1; i < batch_size; ++i) {
                const blst::Fr inv = fr_inverse(fr_from_u64(i));
                kernel[i] = inv;

                blst::Fr neg_inv{};
                blst::Fr::neg(neg_inv, inv);
                kernel[fft_size - i] = neg_inv;
            }
            return kernel;
        }

        std::vector<blst::Fr> build_self_coefficients(const std::size_t batch_size)
        {
            std::vector<blst::Fr> inverse_values(batch_size + 2, fr_zero());
            for (std::size_t i = 1; i <= batch_size + 1; ++i) {
                inverse_values[i] = fr_inverse(fr_from_u64(i));
            }

            std::vector<blst::Fr> harmonic_prefix(batch_size + 1, fr_zero());
            for (std::size_t i = 1; i <= batch_size; ++i) {
                harmonic_prefix[i] = fr_add(harmonic_prefix[i - 1], inverse_values[i]);
            }

            std::vector<blst::Fr> out(batch_size, fr_zero());
            for (std::size_t slot = 1; slot <= batch_size; ++slot) {
                const blst::Fr left = harmonic_prefix[slot - 1];
                const blst::Fr right = harmonic_prefix[batch_size - slot];
                out[slot - 1] = fr_add(
                    fr_add(inverse_values[slot], inverse_values[slot + 1]),
                    fr_sub(left, right));
            }
            return out;
        }

        const PublicParameters& resolve_public_parameters(
            const PreparedPublicParameters& prepared_public_parameters)
        {
            if (prepared_public_parameters.public_parameters == nullptr) {
                throw std::invalid_argument("PFE prepared public parameters are uninitialized");
            }
            return *prepared_public_parameters.public_parameters;
        }

        template<typename G, typename CoeffCtx>
        std::vector<G> toeplitz_multiply_fft_owned(
            const PublicParameters& public_parameters,
            const std::span<const G> input,
            const G& zero)
        {
            if (input.size() != public_parameters.batch_size) {
                throw std::invalid_argument("PFE Toeplitz input has wrong size");
            }

            std::vector<G> fft_buffer(public_parameters.toeplitz_fft_size, zero);
            std::copy(input.begin(), input.end(), fft_buffer.begin());

            fft_forward_truncated_prefix<G, blst::Fr, CoeffCtx>(
                std::span<G>(fft_buffer.data(), fft_buffer.size()),
                std::span<const blst::Fr>(
                    public_parameters.toeplitz_roots.data(),
                    public_parameters.toeplitz_roots.size()),
                input.size());

            for (std::size_t i = 0; i < fft_buffer.size(); ++i) {
                fft_buffer[i] = CoeffCtx::mult(fft_buffer[i], public_parameters.toeplitz_kernel_fft[i]);
            }

            fft_inverse_truncated_range<G, blst::Fr, CoeffCtx>(
                std::span<G>(fft_buffer.data(), fft_buffer.size()),
                std::span<const blst::Fr>(
                    public_parameters.toeplitz_inv_roots.data(),
                    public_parameters.toeplitz_inv_roots.size()),
                0,
                input.size());

            return std::vector<G>(fft_buffer.begin(), fft_buffer.begin() + input.size());
        }

        template<typename G, typename CoeffCtx>
        std::vector<G> toeplitz_multiply_fft_reuse(
            const PublicParameters& public_parameters,
            const std::span<const G> input,
            std::vector<G>& workspace,
            const G& zero)
        {
            if (input.size() != public_parameters.batch_size) {
                throw std::invalid_argument("PFE Toeplitz input has wrong size");
            }

            workspace.assign(public_parameters.toeplitz_fft_size, zero);
            std::copy(input.begin(), input.end(), workspace.begin());

            fft_forward_truncated_prefix<G, blst::Fr, CoeffCtx>(
                std::span<G>(workspace.data(), workspace.size()),
                std::span<const blst::Fr>(
                    public_parameters.toeplitz_roots.data(),
                    public_parameters.toeplitz_roots.size()),
                input.size());

            for (std::size_t i = 0; i < workspace.size(); ++i) {
                workspace[i] = CoeffCtx::mult(workspace[i], public_parameters.toeplitz_kernel_fft[i]);
            }

            fft_inverse_truncated_range<G, blst::Fr, CoeffCtx>(
                std::span<G>(workspace.data(), workspace.size()),
                std::span<const blst::Fr>(
                    public_parameters.toeplitz_inv_roots.data(),
                    public_parameters.toeplitz_inv_roots.size()),
                0,
                input.size());

            return std::vector<G>(workspace.begin(), workspace.begin() + input.size());
        }
    }

    SetupResult Setup(const PfeParameters& params, const RandomBytesSource& rng)
    {
        validate_params(params);

        const blst::Fr x = sample_valid_setup_secret(params, rng);
        PublicParameters public_parameters{};
        public_parameters.batch_size = params.batch_size;
        public_parameters.num_servers = params.num_servers;
        public_parameters.threshold = params.threshold;

        blst::G1Point sum_p = blst::G1Point::infinity();
        for (std::int64_t i = -1; i <= static_cast<std::int64_t>(params.batch_size); ++i) {
            sum_p = add_g1(
                sum_p,
                scale_g1(blst::G1Point::generator(), partial_fraction_value(x, i)));
        }
        public_parameters.sum_p_g1 = sum_p;
        public_parameters.p0_gt = scale_gt(
            pairing(blst::G1Point::generator(), blst::G2Point::generator()),
            partial_fraction_value(x, 0));

        public_parameters.dk.reserve(params.batch_size);
        for (std::size_t slot = 1; slot <= params.batch_size; ++slot) {
            public_parameters.dk.push_back(scale_g2(
                blst::G2Point::generator(),
                partial_fraction_value(x, static_cast<std::int64_t>(slot))));
        }

        public_parameters.sum_dk = blst::G2Point::infinity();
        for (const blst::G2Point& dk_slot : public_parameters.dk) {
            public_parameters.sum_dk = add_g2(public_parameters.sum_dk, dk_slot);
        }

        public_parameters.toeplitz_fft_size = next_power_of_two(params.batch_size * 2);
        public_parameters.toeplitz_roots = make_fft_roots(public_parameters.toeplitz_fft_size);
        public_parameters.toeplitz_inv_roots = make_ifft_roots(std::span<const blst::Fr>(
            public_parameters.toeplitz_roots.data(),
            public_parameters.toeplitz_roots.size()));
        public_parameters.toeplitz_kernel_fft = build_toeplitz_kernel_embedding(
            params.batch_size,
            public_parameters.toeplitz_fft_size);
        fft<blst::Fr, blst::Fr, blst::FrCoeffCtx>(
            std::span<blst::Fr>(
                public_parameters.toeplitz_kernel_fft.data(),
                public_parameters.toeplitz_kernel_fft.size()),
            std::span<const blst::Fr>(
                public_parameters.toeplitz_roots.data(),
                public_parameters.toeplitz_roots.size()),
            false);

        public_parameters.self_coefficients = build_self_coefficients(params.batch_size);
        public_parameters.weighted_dk_excluding_slot = toeplitz_multiply_fft_owned<blst::G2Point, blst::G2CoeffCtx>(
            public_parameters,
            std::span<const blst::G2Point>(public_parameters.dk.data(), public_parameters.dk.size()),
            blst::G2Point::infinity());

        std::vector<ServerSecretShare> secret_shares(params.num_servers);
        for (std::size_t server = 0; server < params.num_servers; ++server) {
            secret_shares[server].index = static_cast<std::uint32_t>(server + 1);
            secret_shares[server].p_shares.assign(params.batch_size, fr_zero());
        }

        for (std::size_t slot = 1; slot <= params.batch_size; ++slot) {
            std::vector<blst::Fr> polynomial(params.threshold, fr_zero());
            polynomial[0] = partial_fraction_value(x, static_cast<std::int64_t>(slot));
            for (std::size_t degree = 1; degree < params.threshold; ++degree) {
                polynomial[degree] = random_nonzero_fr(rng, "BTE_PFE_SETUP_SHARE_V1");
            }

            for (std::size_t server = 0; server < params.num_servers; ++server) {
                secret_shares[server].p_shares[slot - 1] = evaluate_polynomial_at(
                    std::span<const blst::Fr>(polynomial.data(), polynomial.size()),
                    fr_from_u64(server + 1));
            }
        }

        return SetupResult{ std::move(public_parameters), std::move(secret_shares) };
    }

    PreparedPublicParameters PreparePublicParameters(const PublicParameters& public_parameters)
    {
        PreparedPublicParameters prepared{};
        prepared.public_parameters = &public_parameters;
        prepared.g1_fft_workspace.reserve(public_parameters.toeplitz_fft_size);
        prepared.gt_fft_workspace.reserve(public_parameters.toeplitz_fft_size);
        return prepared;
    }

    EncryptionResult Encrypt(
        const PfeParameters& params,
        const PublicParameters& public_parameters,
        const RandomBytesSource& rng)
    {
        validate_params(params);
        if (public_parameters.batch_size != params.batch_size) {
            throw std::invalid_argument("PFE parameters and public parameters disagree on batch_size");
        }

        const blst::Fr r = random_nonzero_fr(rng, "BTE_PFE_ENCRYPT_R_V1");

        EncryptionResult out{};
        out.ciphertext.c1 = scale_g1(blst::G1Point::generator(), r);
        out.ciphertext.c2 = scale_g1(public_parameters.sum_p_g1, r);
        out.key = kdf_from_gt(params, scale_gt(public_parameters.p0_gt, r));
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
            throw std::invalid_argument("PFE server share index must be nonzero");
        }
        if (request.ordered_c1.empty()) {
            throw std::invalid_argument("PFE threshold request must be nonempty");
        }
        if (secret_share.p_shares.size() != request.ordered_c1.size()) {
            throw std::invalid_argument("PFE server share does not match the fixed batch size");
        }

        PartialDecryptionShare share{};
        share.index = secret_share.index;
        share.value = msm_sum_g1(
            std::span<const blst::G1Point>(request.ordered_c1.data(), request.ordered_c1.size()),
            std::span<const blst::Fr>(secret_share.p_shares.data(), secret_share.p_shares.size()));
        return share;
    }

    blst::G1Point CombinePartialDecryptionShares(
        const std::span<const PartialDecryptionShare> shares,
        const PublicParameters& public_parameters)
    {
        if (shares.empty()) {
            throw std::invalid_argument("Need at least one PFE share to combine");
        }
        if (shares.size() < public_parameters.threshold) {
            throw std::invalid_argument("Need at least threshold-many PFE shares to combine");
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
        out.tm_self.resize(ciphertexts.size(), blst::GT::identity());

        for (std::size_t j = 0; j < ciphertexts.size(); ++j) {
            out.tm_self[j] = pairing(ciphertexts[j].c1, public_parameters.dk[j]);
        }

        std::vector<blst::G1Point> c1_points;
        c1_points.reserve(ciphertexts.size());
        for (const Ciphertext& ciphertext : ciphertexts) {
            c1_points.push_back(ciphertext.c1);
        }

        out.weighted_c1_excluding_slot = toeplitz_multiply_fft_owned<blst::G1Point, blst::G1CoeffCtx>(
            public_parameters,
            std::span<const blst::G1Point>(c1_points.data(), c1_points.size()),
            blst::G1Point::infinity());
        out.weighted_tm_excluding_slot = toeplitz_multiply_fft_owned<blst::GT, blst::GTCoeffCtx>(
            public_parameters,
            std::span<const blst::GT>(out.tm_self.data(), out.tm_self.size()),
            blst::GT::identity());

        return out;
    }

    BatchDecryptionPrecomputation PrecomputeBatchDecryption(
        const PreparedPublicParameters& prepared_public_parameters,
        const std::span<const Ciphertext> ciphertexts)
    {
        const PublicParameters& public_parameters = resolve_public_parameters(prepared_public_parameters);
        validate_batch_size(public_parameters, ciphertexts.size());

        BatchDecryptionPrecomputation out{};
        out.ciphertexts.assign(ciphertexts.begin(), ciphertexts.end());
        out.tm_self.resize(ciphertexts.size(), blst::GT::identity());

        for (std::size_t j = 0; j < ciphertexts.size(); ++j) {
            out.tm_self[j] = pairing(ciphertexts[j].c1, public_parameters.dk[j]);
        }

        std::vector<blst::G1Point> c1_points;
        c1_points.reserve(ciphertexts.size());
        for (const Ciphertext& ciphertext : ciphertexts) {
            c1_points.push_back(ciphertext.c1);
        }

        out.weighted_c1_excluding_slot = toeplitz_multiply_fft_reuse<blst::G1Point, blst::G1CoeffCtx>(
            public_parameters,
            std::span<const blst::G1Point>(c1_points.data(), c1_points.size()),
            prepared_public_parameters.g1_fft_workspace,
            blst::G1Point::infinity());
        out.weighted_tm_excluding_slot = toeplitz_multiply_fft_reuse<blst::GT, blst::GTCoeffCtx>(
            public_parameters,
            std::span<const blst::GT>(out.tm_self.data(), out.tm_self.size()),
            prepared_public_parameters.gt_fft_workspace,
            blst::GT::identity());

        return out;
    }

    std::vector<Bytes> FinalizeBatchDecryption(
        const PfeParameters& params,
        const PublicParameters& public_parameters,
        const BatchDecryptionPrecomputation& precomputation,
        const blst::G1Point& aggregate_predecryption_key)
    {
        validate_params(params);
        validate_batch_size(public_parameters, precomputation.ciphertexts.size());

        if (precomputation.tm_self.size() != public_parameters.batch_size ||
            precomputation.weighted_c1_excluding_slot.size() != public_parameters.batch_size ||
            precomputation.weighted_tm_excluding_slot.size() != public_parameters.batch_size) {
            throw std::invalid_argument("PFE precomputation vectors have the wrong size");
        }

        std::vector<Bytes> keys;
        keys.reserve(public_parameters.batch_size);
        for (std::size_t j = 0; j < public_parameters.batch_size; ++j) {
            const Ciphertext& ciphertext = precomputation.ciphertexts[j];

            blst::G1Point l1_input = sub_g1(ciphertext.c2, aggregate_predecryption_key);
            l1_input = add_g1(
                l1_input,
                scale_g1(ciphertext.c1, public_parameters.self_coefficients[j]));
            l1_input = sub_g1(l1_input, precomputation.weighted_c1_excluding_slot[j]);

            blst::GT l1 = pairing(l1_input, public_parameters.dk[j]);
            l1 = sub_gt(
                l1,
                pairing(ciphertext.c1, public_parameters.weighted_dk_excluding_slot[j]));
            l1 = add_gt(l1, precomputation.weighted_tm_excluding_slot[j]);

            blst::GT l2 = pairing(ciphertext.c2, blst::G2Point::generator());
            l2 = sub_gt(l2, pairing(ciphertext.c1, public_parameters.sum_dk));

            const blst::Fr slot = fr_from_u64(j + 1);
            const blst::Fr slot_plus_one = fr_from_u64(j + 2);
            const blst::GT rp0 = scale_gt(
                sub_gt(scale_gt(l1, slot_plus_one), l2),
                slot);
            keys.push_back(kdf_from_gt(params, rp0));
        }

        return keys;
    }

    std::vector<Bytes> FinalizeBatchDecryption(
        const PfeParameters& params,
        const PreparedPublicParameters& prepared_public_parameters,
        const BatchDecryptionPrecomputation& precomputation,
        const blst::G1Point& aggregate_predecryption_key)
    {
        return FinalizeBatchDecryption(
            params,
            resolve_public_parameters(prepared_public_parameters),
            precomputation,
            aggregate_predecryption_key);
    }

    std::vector<Bytes> RecoverBatchKeys(
        const PfeParameters& params,
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
        const PfeParameters& params,
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
