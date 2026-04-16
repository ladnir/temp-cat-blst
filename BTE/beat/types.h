#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <span>
#include <string>
#include <vector>

#include "blst/curve/GT.h"
#include "blst/curve/point.h"
#include "blst/fields/Fr.h"
#include "blst/fft/middle_product.h"

namespace bte::beat
{
    using Bytes = std::vector<std::uint8_t>;
    using RandomBytesSource = std::function<void(std::span<std::uint8_t>)>;

    struct BeatParameters
    {
        std::size_t max_batch_size = 0;
        std::size_t num_servers = 0;
        std::size_t threshold = 0;
        std::size_t key_size = 32;

        std::string kdf_salt = "BTE_BEAT_KDF_SALT_V1";
        std::string kdf_info = "BTE_BEAT_KDF_INFO_V1";
    };

    struct PublicParameters
    {
        std::size_t max_batch_size = 0;
        std::size_t num_servers = 0;
        std::size_t threshold = 0;
        blst::GT gt_tau_to_bmax_plus_one = blst::GT::identity();
        std::vector<blst::G2Point> g2_powers;
    };

    struct ServerSecretShare
    {
        std::uint32_t index = 0;
        std::vector<blst::Fr> tau_power_shares;
    };

    struct SetupResult
    {
        PublicParameters public_parameters;
        std::vector<ServerSecretShare> secret_shares;
    };

    struct PreparedPublicParameters
    {
        PublicParameters public_parameters;
        mutable std::vector<blst::MiddleProductPrecomputedB> beta_precomputed_b_by_batch_size;
        mutable std::vector<blst::MiddleProductPrecomputedB> alpha_precomputed_b_by_batch_size;
    };

    struct Ciphertext
    {
        blst::G1Point c1 = blst::G1Point::infinity();
    };

    struct EncryptionResult
    {
        Bytes key;
        Ciphertext ciphertext;
    };

    struct ThresholdRequest
    {
        std::vector<blst::G1Point> ordered_c1;
    };

    struct PartialDecryptionShare
    {
        std::uint32_t index = 0;
        blst::G1Point value = blst::G1Point::infinity();
    };

    struct BatchDecryptionPrecomputation
    {
        std::vector<Ciphertext> ciphertexts;
        std::vector<blst::GT> beta;
    };
}
