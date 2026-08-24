#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "blst/curve/point.h"
#include "blst/fields/Fr.h"

namespace bte::btx
{
    enum class CiphertextAdmissionError
    {
        Ok,
        InvalidLength,
        InvalidPoint,
        InvalidScalar,
        IdentityPoint,
        ClientNizkFailed,
    };

    struct ClientProof
    {
        blst::Fr challenge = blst::Fr::zero();
        blst::Fr response = blst::Fr::zero();
    };

    // This is the protocol ciphertext admitted from the transaction wire. The
    // legacy Ciphertext in types.h remains the KEM point used by the existing
    // C++ decryption prototype until the two internal representations are unified.
    struct ClientCiphertext
    {
        blst::G1Point r = blst::G1Point::infinity();
        std::array<std::uint8_t, 16> c1{};
        std::vector<std::uint8_t> c2;
        ClientProof proof{};
    };

    // Decodes R (48) || C1 (16) || len(C2) (4-byte big endian) || C2 ||
    // challenge (32-byte canonical big endian) || response (same).
    // Point decoding establishes canonical on-curve encoding only. Verification
    // performs the subgroup and identity checks before any scalar multiplication.
    CiphertextAdmissionError DecodeClientCiphertext(
        std::span<const std::uint8_t> encoded,
        ClientCiphertext& out);

    CiphertextAdmissionError VerifyClientCiphertext(
        const ClientCiphertext& ciphertext,
        std::span<const std::uint8_t> associated_data);

    CiphertextAdmissionError DecodeAndVerifyClientCiphertext(
        std::span<const std::uint8_t> encoded,
        std::span<const std::uint8_t> associated_data,
        ClientCiphertext& out);

    const char* CiphertextAdmissionErrorId(CiphertextAdmissionError error);
}
