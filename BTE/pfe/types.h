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

namespace bte::pfe
{
    using Bytes = std::vector<std::uint8_t>;
    using RandomBytesSource = std::function<void(std::span<std::uint8_t>)>;

    struct PfeParameters
    {
        std::size_t batch_size = 0;
        std::size_t num_servers = 0;
        std::size_t threshold = 0;
        std::size_t key_size = 32;

        std::string kdf_salt = "BTE_PFE_KDF_SALT_V1";
        std::string kdf_info = "BTE_PFE_KDF_INFO_V1";
    };

    struct PublicParameters
    {
        std::size_t batch_size = 0;
        std::size_t num_servers = 0;
        std::size_t threshold = 0;

        blst::G1Point sum_p_g1 = blst::G1Point::infinity();
        blst::GT p0_gt = blst::GT::identity();

        std::vector<blst::G2Point> dk;
        blst::G2Point sum_dk = blst::G2Point::infinity();
        std::vector<blst::G2Point> weighted_dk_excluding_slot;
        std::vector<blst::Fr> self_coefficients;
        std::size_t toeplitz_fft_size = 0;
        std::vector<blst::Fr> toeplitz_roots;
        std::vector<blst::Fr> toeplitz_inv_roots;
        std::vector<blst::Fr> toeplitz_kernel_fft;
    };

    struct ServerSecretShare
    {
        std::uint32_t index = 0;
        std::vector<blst::Fr> p_shares;
    };

    struct SetupResult
    {
        PublicParameters public_parameters;
        std::vector<ServerSecretShare> secret_shares;
    };

    struct PreparedPublicParameters
    {
        const PublicParameters* public_parameters = nullptr;
        mutable std::vector<blst::G1Point> g1_fft_workspace;
        mutable std::vector<blst::GT> gt_fft_workspace;
    };

    struct Ciphertext
    {
        blst::G1Point c1 = blst::G1Point::infinity();
        blst::G1Point c2 = blst::G1Point::infinity();
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
        std::vector<blst::GT> tm_self;
        std::vector<blst::G1Point> weighted_c1_excluding_slot;
        std::vector<blst::GT> weighted_tm_excluding_slot;
    };
}
