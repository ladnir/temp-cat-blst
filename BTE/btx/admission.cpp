#include "BTE/btx/admission.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>

#include "bindings/blst.h"
#include "blake3.h"

namespace bte::btx
{
    namespace
    {
        constexpr std::size_t kPointSize = 48;
        constexpr std::size_t kC1Size = 16;
        constexpr std::size_t kLengthSize = 4;
        constexpr std::size_t kScalarSize = 32;
        constexpr std::size_t kProofSize = 2 * kScalarSize;
        constexpr std::size_t kMinimumCiphertextSize =
            kPointSize + kC1Size + kLengthSize + kProofSize;

        std::uint32_t read_u32_be(const std::uint8_t* bytes)
        {
            return (static_cast<std::uint32_t>(bytes[0]) << 24) |
                (static_cast<std::uint32_t>(bytes[1]) << 16) |
                (static_cast<std::uint32_t>(bytes[2]) << 8) |
                static_cast<std::uint32_t>(bytes[3]);
        }

        bool decode_scalar_be(const std::uint8_t* bytes, blst::Fr& out)
        {
            blst_scalar scalar{};
            blst_scalar_from_bendian(&scalar, bytes);
            if (!blst_scalar_fr_check(&scalar)) {
                return false;
            }
            blst_fr_from_scalar(&out, &scalar);
            return true;
        }

        void update_length_prefixed(
            blake3_hasher& hasher,
            const std::span<const std::uint8_t> bytes)
        {
            std::array<std::uint8_t, 8> length{};
            std::uint64_t value = static_cast<std::uint64_t>(bytes.size());
            for (std::size_t i = 0; i < length.size(); ++i) {
                length[length.size() - 1 - i] = static_cast<std::uint8_t>(value);
                value >>= 8;
            }
            blake3_hasher_update(&hasher, length.data(), length.size());
            blake3_hasher_update(&hasher, bytes.data(), bytes.size());
        }

        blst::Fr challenge(
            const blst::G1Point& r,
            const blst::G1Point& t,
            const std::array<std::uint8_t, kC1Size>& c1,
            const std::span<const std::uint8_t> c2,
            const std::span<const std::uint8_t> associated_data)
        {
            std::array<std::uint8_t, kPointSize> compressed_r{};
            std::array<std::uint8_t, kPointSize> compressed_t{};
            blst::G1Point::compress(compressed_r.data(), r);
            blst::G1Point::compress(compressed_t.data(), t);

            blake3_hasher hasher{};
            blake3_hasher_init_derive_key(&hasher, "btx/nizk/v1");
            blake3_hasher_update(&hasher, compressed_r.data(), compressed_r.size());
            blake3_hasher_update(&hasher, compressed_t.data(), compressed_t.size());
            blake3_hasher_update(&hasher, c1.data(), c1.size());
            update_length_prefixed(hasher, c2);
            update_length_prefixed(hasher, associated_data);

            std::array<std::uint8_t, 64> wide{};
            blake3_hasher_finalize(&hasher, wide.data(), wide.size());
            blst_scalar scalar{};
            (void)blst_scalar_from_be_bytes(&scalar, wide.data(), wide.size());
            blst::Fr result{};
            blst_fr_from_scalar(&result, &scalar);
            return result;
        }

        void public_scalar_multiply(
            blst::G1Point& out,
            const blst::G1Point& point,
            const blst::Fr& scalar)
        {
            std::array<std::uint8_t, kScalarSize> scalar_le{};
            blst::Fr::to_scaler(scalar_le.data(), scalar);
            blst::G1Point::mult_glsv_vartime(out, point, scalar_le.data());
        }
    }

    CiphertextAdmissionError DecodeClientCiphertext(
        const std::span<const std::uint8_t> encoded,
        ClientCiphertext& out)
    {
        if (encoded.size() < kMinimumCiphertextSize) {
            return CiphertextAdmissionError::InvalidLength;
        }

        blst::G1AffinePoint affine{};
        if (blst::G1AffinePoint::uncompress(affine, encoded.data()) != blst::BLST_SUCCESS) {
            return CiphertextAdmissionError::InvalidPoint;
        }

        const std::size_t c2_size = read_u32_be(encoded.data() + kPointSize + kC1Size);
        if (c2_size > std::numeric_limits<std::size_t>::max() - kMinimumCiphertextSize ||
            encoded.size() != kMinimumCiphertextSize + c2_size) {
            return CiphertextAdmissionError::InvalidLength;
        }

        ClientCiphertext decoded{};
        blst::G1Point::from_affine(decoded.r, affine);
        std::copy_n(encoded.data() + kPointSize, kC1Size, decoded.c1.begin());
        const std::uint8_t* const c2_begin =
            encoded.data() + kPointSize + kC1Size + kLengthSize;
        decoded.c2.assign(c2_begin, c2_begin + c2_size);
        const std::uint8_t* const proof = c2_begin + c2_size;
        if (!decode_scalar_be(proof, decoded.proof.challenge) ||
            !decode_scalar_be(proof + kScalarSize, decoded.proof.response)) {
            return CiphertextAdmissionError::InvalidScalar;
        }

        out = std::move(decoded);
        return CiphertextAdmissionError::Ok;
    }

    CiphertextAdmissionError VerifyClientCiphertext(
        const ClientCiphertext& ciphertext,
        const std::span<const std::uint8_t> associated_data)
    {
        if (blst::G1Point::in_subgroup(ciphertext.r) == 0) {
            return CiphertextAdmissionError::InvalidPoint;
        }
        if (blst::G1Point::is_identity(ciphertext.r) != 0) {
            return CiphertextAdmissionError::IdentityPoint;
        }

        blst::G1Point s_g{};
        blst::G1Point c_r{};
        public_scalar_multiply(s_g, blst::G1Point::generator(), ciphertext.proof.response);
        public_scalar_multiply(c_r, ciphertext.r, ciphertext.proof.challenge);
        blst::G1Point t{};
        blst::G1Point::dadd(t, s_g, c_r);

        const blst::Fr expected = challenge(
            ciphertext.r,
            t,
            ciphertext.c1,
            ciphertext.c2,
            associated_data);
        return expected == ciphertext.proof.challenge
            ? CiphertextAdmissionError::Ok
            : CiphertextAdmissionError::ClientNizkFailed;
    }

    CiphertextAdmissionError DecodeAndVerifyClientCiphertext(
        const std::span<const std::uint8_t> encoded,
        const std::span<const std::uint8_t> associated_data,
        ClientCiphertext& out)
    {
        ClientCiphertext decoded{};
        const CiphertextAdmissionError decode_error = DecodeClientCiphertext(encoded, decoded);
        if (decode_error != CiphertextAdmissionError::Ok) {
            return decode_error;
        }
        const CiphertextAdmissionError verify_error =
            VerifyClientCiphertext(decoded, associated_data);
        if (verify_error != CiphertextAdmissionError::Ok) {
            return verify_error;
        }
        out = std::move(decoded);
        return CiphertextAdmissionError::Ok;
    }

    const char* CiphertextAdmissionErrorId(const CiphertextAdmissionError error)
    {
        switch (error) {
        case CiphertextAdmissionError::Ok:
            return "ok";
        case CiphertextAdmissionError::InvalidLength:
            return "invalid_length";
        case CiphertextAdmissionError::InvalidPoint:
            return "invalid_point";
        case CiphertextAdmissionError::InvalidScalar:
            return "invalid_scalar";
        case CiphertextAdmissionError::IdentityPoint:
            return "identity_point";
        case CiphertextAdmissionError::ClientNizkFailed:
            return "client_nizk_failed";
        }
        return "unknown_admission_error";
    }
}
