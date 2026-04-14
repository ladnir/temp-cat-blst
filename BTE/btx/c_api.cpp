#include "bindings/blst_btx.h"

#include <cstring>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "BTE/btx/core.h"
#include "blst/scalar/curve/point.h"
#if BLST_DEFAULT_AVX
#include "blst/avx/curve/map_to_G.h"
#endif

namespace
{
    using BoundaryG1 = blst::scalar::G1Point;

    static_assert(sizeof(blst_p1) == sizeof(BoundaryG1));
    static_assert(alignof(blst_p1) == alignof(BoundaryG1));

    bte::btx::BtxParameters to_cpp_params(const blst_btx_params& params)
    {
        bte::btx::BtxParameters out{};
        out.max_batch_size = params.max_batch_size;
        out.num_servers = params.num_servers;
        out.threshold = params.threshold;
        out.key_size = params.key_size;
        if (params.kdf_salt != nullptr) {
            out.kdf_salt = params.kdf_salt;
        }
        if (params.kdf_info != nullptr) {
            out.kdf_info = params.kdf_info;
        }
        return out;
    }

    bte::btx::RandomBytesSource to_rng(
        blst_fill_random_fn fill_random,
        void* user_data)
    {
        return [fill_random, user_data](std::span<std::uint8_t> out) {
            if (fill_random == nullptr || !fill_random(user_data, out.data(), out.size())) {
                throw std::runtime_error("BTX random callback failed");
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

    bte::btx::CrossTermAlgorithm to_cpp_algorithm(
        const blst_btx_cross_term_algorithm algorithm)
    {
        switch (algorithm) {
        case BLST_BTX_CROSS_TERM_NAIVE:
            return bte::btx::CrossTermAlgorithm::Naive;
        case BLST_BTX_CROSS_TERM_FFT:
            return bte::btx::CrossTermAlgorithm::MiddleProductFft;
        case BLST_BTX_CROSS_TERM_MULTIPAIR:
            return bte::btx::CrossTermAlgorithm::MiddleProductMultipair;
        case BLST_BTX_CROSS_TERM_AUTO:
            return bte::btx::CrossTermAlgorithm::MiddleProductAuto;
        default:
            throw std::invalid_argument("Unknown BTX cross-term algorithm");
        }
    }

    std::vector<bte::btx::Ciphertext> to_ciphertexts(
        const blst_p1* ciphertexts,
        const std::size_t n)
    {
        std::vector<bte::btx::Ciphertext> out;
        out.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            out.push_back({ to_public_g1(ciphertexts[i]) });
        }
        return out;
    }

    bool copy_keys(
        const std::vector<bte::btx::Bytes>& keys,
        const std::size_t key_size,
        byte* out_keys,
        const std::size_t out_keys_len)
    {
        if (out_keys == nullptr || out_keys_len != (keys.size() * key_size)) {
            return false;
        }
        for (std::size_t i = 0; i < keys.size(); ++i) {
            if (keys[i].size() != key_size) {
                return false;
            }
            std::memcpy(out_keys + (i * key_size), keys[i].data(), key_size);
        }
        return true;
    }
}

struct blst_btx_public_parameters
{
    const bte::btx::PublicParameters* value = nullptr;
};

struct blst_btx_server_secret_share
{
    const bte::btx::ServerSecretShare* value = nullptr;
};

struct blst_btx_setup_result
{
    bte::btx::SetupResult value;
    blst_btx_public_parameters public_parameters_view{};
    std::vector<blst_btx_server_secret_share> secret_share_views;
};

struct blst_btx_prepared_public_parameters
{
    bte::btx::PreparedPublicParameters value;
};

struct blst_btx_batch_precomputation
{
    bte::btx::BatchDecryptionPrecomputation value;
};

extern "C" {

blst_btx_setup_result* blst_btx_setup(
    const blst_btx_params* params,
    blst_fill_random_fn fill_random,
    void* user_data)
{
    if (params == nullptr || fill_random == nullptr) {
        return nullptr;
    }

    try {
        auto out = std::make_unique<blst_btx_setup_result>();
        out->value = bte::btx::Setup(to_cpp_params(*params), to_rng(fill_random, user_data));
        out->public_parameters_view.value = &out->value.public_parameters;
        out->secret_share_views.resize(out->value.secret_shares.size());
        for (std::size_t i = 0; i < out->value.secret_shares.size(); ++i) {
            out->secret_share_views[i].value = &out->value.secret_shares[i];
        }
        return out.release();
    } catch (...) {
        return nullptr;
    }
}

const blst_btx_public_parameters* blst_btx_setup_result_public_parameters(
    const blst_btx_setup_result* setup_result)
{
    if (setup_result == nullptr) {
        return nullptr;
    }
    return &setup_result->public_parameters_view;
}

size_t blst_btx_setup_result_secret_share_count(
    const blst_btx_setup_result* setup_result)
{
    if (setup_result == nullptr) {
        return 0;
    }
    return setup_result->secret_share_views.size();
}

const blst_btx_server_secret_share* blst_btx_setup_result_secret_share(
    const blst_btx_setup_result* setup_result,
    const size_t index)
{
    if (setup_result == nullptr || index >= setup_result->secret_share_views.size()) {
        return nullptr;
    }
    return &setup_result->secret_share_views[index];
}

void blst_btx_setup_result_free(blst_btx_setup_result* setup_result)
{
    delete setup_result;
}

blst_btx_prepared_public_parameters* blst_btx_prepare_public_parameters(
    const blst_btx_public_parameters* public_parameters)
{
    if (public_parameters == nullptr || public_parameters->value == nullptr) {
        return nullptr;
    }

    try {
        auto out = std::make_unique<blst_btx_prepared_public_parameters>();
        out->value = bte::btx::PreparePublicParameters(*public_parameters->value);
        return out.release();
    } catch (...) {
        return nullptr;
    }
}

void blst_btx_prepared_public_parameters_free(
    blst_btx_prepared_public_parameters* prepared_public_parameters)
{
    delete prepared_public_parameters;
}

bool blst_btx_encrypt(
    const blst_btx_params* params,
    const blst_btx_public_parameters* public_parameters,
    blst_fill_random_fn fill_random,
    void* user_data,
    blst_p1* out_c1,
    byte* out_key,
    const size_t out_key_len)
{
    if (params == nullptr || public_parameters == nullptr || public_parameters->value == nullptr ||
        fill_random == nullptr || out_c1 == nullptr || out_key == nullptr) {
        return false;
    }

    try {
        const bte::btx::EncryptionResult encrypted = bte::btx::Encrypt(
            to_cpp_params(*params),
            *public_parameters->value,
            to_rng(fill_random, user_data));
        if (encrypted.key.size() != out_key_len) {
            return false;
        }
        *out_c1 = to_boundary_g1(encrypted.ciphertext.c1);
        std::memcpy(out_key, encrypted.key.data(), encrypted.key.size());
        return true;
    } catch (...) {
        return false;
    }
}

bool blst_btx_create_partial_decryption_share(
    const blst_btx_server_secret_share* secret_share,
    const blst_p1* ordered_c1,
    const size_t n,
    blst_btx_partial_decryption_share* out_share)
{
    if (secret_share == nullptr || secret_share->value == nullptr ||
        ordered_c1 == nullptr || n == 0 || out_share == nullptr) {
        return false;
    }

    try {
        std::vector<blst::G1Point> ordered;
        ordered.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            ordered.push_back(to_public_g1(ordered_c1[i]));
        }
        const bte::btx::ThresholdRequest request{ std::move(ordered) };
        const bte::btx::PartialDecryptionShare share =
            bte::btx::CreatePartialDecryptionShare(*secret_share->value, request);
        out_share->index = share.index;
        out_share->value = to_boundary_g1(share.value);
        return true;
    } catch (...) {
        return false;
    }
}

bool blst_btx_combine_partial_decryption_shares(
    const blst_btx_partial_decryption_share* shares,
    const size_t n,
    const blst_btx_public_parameters* public_parameters,
    blst_p1* out_aggregate_predecryption_key)
{
    if (shares == nullptr || n == 0 || public_parameters == nullptr ||
        public_parameters->value == nullptr || out_aggregate_predecryption_key == nullptr) {
        return false;
    }

    try {
        std::vector<bte::btx::PartialDecryptionShare> cpp_shares;
        cpp_shares.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            cpp_shares.push_back({ shares[i].index, to_public_g1(shares[i].value) });
        }

        const blst::G1Point aggregate = bte::btx::CombinePartialDecryptionShares(
            std::span<const bte::btx::PartialDecryptionShare>(cpp_shares.data(), cpp_shares.size()),
            *public_parameters->value);
        *out_aggregate_predecryption_key = to_boundary_g1(aggregate);
        return true;
    } catch (...) {
        return false;
    }
}

blst_btx_batch_precomputation* blst_btx_precompute_batch_decryption(
    const blst_btx_prepared_public_parameters* prepared_public_parameters,
    const blst_p1* ciphertexts,
    const size_t n,
    const blst_btx_cross_term_algorithm algorithm)
{
    if (prepared_public_parameters == nullptr || ciphertexts == nullptr || n == 0) {
        return nullptr;
    }

    try {
        const std::vector<bte::btx::Ciphertext> cpp_ciphertexts = to_ciphertexts(ciphertexts, n);
        auto out = std::make_unique<blst_btx_batch_precomputation>();
        out->value = bte::btx::PrecomputeBatchDecryption(
            prepared_public_parameters->value,
            std::span<const bte::btx::Ciphertext>(cpp_ciphertexts.data(), cpp_ciphertexts.size()),
            to_cpp_algorithm(algorithm));
        return out.release();
    } catch (...) {
        return nullptr;
    }
}

size_t blst_btx_batch_precomputation_batch_size(
    const blst_btx_batch_precomputation* batch_precomputation)
{
    if (batch_precomputation == nullptr) {
        return 0;
    }
    return batch_precomputation->value.ciphertexts.size();
}

void blst_btx_batch_precomputation_free(
    blst_btx_batch_precomputation* batch_precomputation)
{
    delete batch_precomputation;
}

bool blst_btx_finalize_batch_decryption(
    const blst_btx_params* params,
    const blst_btx_prepared_public_parameters* prepared_public_parameters,
    const blst_btx_batch_precomputation* batch_precomputation,
    const blst_p1* aggregate_predecryption_key,
    byte* out_keys,
    const size_t out_keys_len)
{
    if (params == nullptr || prepared_public_parameters == nullptr ||
        batch_precomputation == nullptr || aggregate_predecryption_key == nullptr) {
        return false;
    }

    try {
        const std::vector<bte::btx::Bytes> keys = bte::btx::FinalizeBatchDecryption(
            to_cpp_params(*params),
            prepared_public_parameters->value,
            batch_precomputation->value,
            to_public_g1(*aggregate_predecryption_key));
        return copy_keys(keys, params->key_size, out_keys, out_keys_len);
    } catch (...) {
        return false;
    }
}

bool blst_btx_recover_batch_keys(
    const blst_btx_params* params,
    const blst_btx_prepared_public_parameters* prepared_public_parameters,
    const blst_p1* ciphertexts,
    const size_t n,
    const blst_p1* aggregate_predecryption_key,
    const blst_btx_cross_term_algorithm algorithm,
    byte* out_keys,
    const size_t out_keys_len)
{
    if (params == nullptr || prepared_public_parameters == nullptr ||
        ciphertexts == nullptr || n == 0 || aggregate_predecryption_key == nullptr) {
        return false;
    }

    try {
        const std::vector<bte::btx::Ciphertext> cpp_ciphertexts = to_ciphertexts(ciphertexts, n);
        const std::vector<bte::btx::Bytes> keys = bte::btx::RecoverBatchKeys(
            to_cpp_params(*params),
            prepared_public_parameters->value,
            std::span<const bte::btx::Ciphertext>(cpp_ciphertexts.data(), cpp_ciphertexts.size()),
            to_public_g1(*aggregate_predecryption_key),
            to_cpp_algorithm(algorithm));
        return copy_keys(keys, params->key_size, out_keys, out_keys_len);
    } catch (...) {
        return false;
    }
}

}
