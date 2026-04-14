#include "BTE/building_blocks/blst_cpp.h"

#include <stdexcept>
#include <vector>

#include "BTE/building_blocks/serialization.h"
#include "blst/curve/map_to_G.h"
#include "blst/fields/Fr.h"
#include "blst/fields/hash_to_field.h"

namespace bte::building_blocks
{
    namespace
    {
        std::array<std::uint8_t, 32> expand_to_valid_scalar(
            const std::span<const std::uint8_t> message,
            const std::string& dst)
        {
            for (std::uint32_t counter = 0; counter < 256; ++counter) {
                std::vector<std::uint8_t> input(message.begin(), message.end());
                append_u32_be(input, counter);

                std::array<std::uint8_t, 32> candidate{};
                blst::blst_expand_message_xmd(
                    candidate.data(),
                    candidate.size(),
                    input.data(),
                    input.size(),
                    reinterpret_cast<const unsigned char*>(dst.data()),
                    dst.size());

                if (blst::sk_check(candidate.data()) != 0) {
                    return candidate;
                }
            }

            throw std::runtime_error("Unable to derive a valid nonzero scalar");
        }
    }

    std::array<std::uint8_t, 32> hash_to_scalar_le(
        const std::span<const std::uint8_t> message,
        const std::string& dst)
    {
        return expand_to_valid_scalar(message, dst);
    }

    blst::G1Point hash_to_g1(
        const std::span<const std::uint8_t> message,
        const std::string& dst)
    {
        blst::G1Point point{};
        blst::Hash_to_curve<blst::G1Point>(
            point,
            reinterpret_cast<const unsigned char*>(message.data()),
            message.size(),
            reinterpret_cast<const unsigned char*>(dst.data()),
            dst.size(),
            nullptr,
            0);
        return point;
    }

    blst::G2Point hash_to_g2(
        const std::span<const std::uint8_t> message,
        const std::string& dst)
    {
        blst::G2Point point{};
        blst::Hash_to_curve<blst::G2Point>(
            point,
            reinterpret_cast<const unsigned char*>(message.data()),
            message.size(),
            reinterpret_cast<const unsigned char*>(dst.data()),
            dst.size(),
            nullptr,
            0);
        return point;
    }
}
