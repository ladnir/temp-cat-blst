#include "BTE/building_blocks/bls_ots.h"

#include <cstring>
#include <stdexcept>

#include "BTE/building_blocks/blst_cpp.h"
#include "blst/curve/GT.h"
#include "blst/fields/Fr.h"

namespace bte::building_blocks
{
    namespace
    {
        bool fr_is_zero(const blst::Fr& value)
        {
            blst::Fr zero{};
            std::memset(&zero, 0, sizeof(zero));
            return blst::Fr::equal(value, zero);
        }

        blst::Fr scalar_bytes_to_fr(const std::array<std::uint8_t, 32>& bytes)
        {
            blst::Fr out{};
            blst::Fr::from_scalar(out, bytes.data());
            return out;
        }

        std::array<std::uint8_t, 32> fr_to_scalar_bytes(const blst::Fr& scalar)
        {
            std::array<std::uint8_t, 32> bytes{};
            blst::Fr::to_scaler(bytes.data(), scalar);
            return bytes;
        }

        blst::G1Point scale_g1(const blst::G1Point& point, const blst::Fr& scalar)
        {
            const std::array<std::uint8_t, 32> scalar_bytes = fr_to_scalar_bytes(scalar);
            blst::G1Point out{};
            blst::G1Point::mult(out, point, scalar_bytes.data(), 256);
            return out;
        }

        blst::G2Point scale_g2(const blst::G2Point& point, const blst::Fr& scalar)
        {
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

        blst::G1Point deserialize_g1(const std::span<const std::uint8_t> bytes)
        {
            if (bytes.size() != 48) {
                throw std::invalid_argument("Compressed G1 point must be 48 bytes");
            }

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

        blst::G2Point deserialize_g2(const std::span<const std::uint8_t> bytes)
        {
            if (bytes.size() != 96) {
                throw std::invalid_argument("Compressed G2 point must be 96 bytes");
            }

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
    }

    BlsMinPkOneTimeSignature sign_bls_min_pk(
        const std::span<const std::uint8_t> ikm,
        const std::span<const std::uint8_t> message,
        const std::string& dst)
    {
        if (ikm.empty()) {
            throw std::invalid_argument("BLS one-time signing key material must not be empty");
        }

        blst::Fr signing_key = scalar_bytes_to_fr(hash_to_scalar_le(ikm, "BTE_TRX_OTS_SK_V1"));
        if (fr_is_zero(signing_key)) {
            throw std::runtime_error("Derived one-time signing key is zero");
        }

        BlsMinPkOneTimeSignature output{};
        const blst::G1Point verification_key = scale_g1(blst::G1Point::generator(), signing_key);
        blst::G1Point::compress(output.verification_key.data(), verification_key);

        const blst::G2Point hashed_message = hash_to_g2(message, dst);
        const blst::G2Point signature = scale_g2(hashed_message, signing_key);
        blst::G2Point::compress(output.signature.data(), signature);
        return output;
    }

    bool verify_bls_min_pk(
        const std::span<const std::uint8_t> verification_key,
        const std::span<const std::uint8_t> signature,
        const std::span<const std::uint8_t> message,
        const std::string& dst)
    {
        if (verification_key.size() != 48 || signature.size() != 96) {
            return false;
        }

        try {
            const blst::G1Point vk = deserialize_g1(verification_key);
            const blst::G2Point sig = deserialize_g2(signature);
            const blst::G2Point hashed_message = hash_to_g2(message, dst);

            return blst::GT::equal(pairing(vk, hashed_message), pairing(blst::G1Point::generator(), sig)) != 0;
        } catch (...) {
            return false;
        }
    }
}
