#include "BTE/btx/serialization.h"

#include <array>
#include <stdexcept>

#include "BTE/building_blocks/serialization.h"

namespace bte::btx
{
    namespace
    {
        constexpr std::uint8_t kBtxWireVersion = 1;

        using Bytes = bte::building_blocks::Bytes;

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

        std::array<std::uint8_t, 12 * 48> serialize_gt(const blst::GT& value)
        {
            std::array<std::uint8_t, 12 * 48> out{};
            blst::GT::serialize(out.data(), value);
            return out;
        }

        blst::G1Point deserialize_g1(const std::span<const std::uint8_t> bytes)
        {
            if (bytes.size() != 48) {
                throw std::invalid_argument("BTX G1 encoding must be 48 bytes");
            }

            blst::G1AffinePoint affine{};
            if (blst::G1AffinePoint::uncompress(affine, bytes.data()) != blst::BLST_SUCCESS) {
                throw std::invalid_argument("Invalid BTX compressed G1 point");
            }
            if (blst::G1AffinePoint::in_subgroup(affine) == 0) {
                throw std::invalid_argument("BTX compressed G1 point is not in the subgroup");
            }

            blst::G1Point out{};
            blst::G1Point::from_affine(out, affine);
            return out;
        }

        blst::G2Point deserialize_g2(const std::span<const std::uint8_t> bytes)
        {
            if (bytes.size() != 96) {
                throw std::invalid_argument("BTX G2 encoding must be 96 bytes");
            }

            blst::G2AffinePoint affine{};
            if (blst::G2AffinePoint::uncompress(affine, bytes.data()) != blst::BLST_SUCCESS) {
                throw std::invalid_argument("Invalid BTX compressed G2 point");
            }
            if (blst::G2AffinePoint::in_subgroup(affine) == 0) {
                throw std::invalid_argument("BTX compressed G2 point is not in the subgroup");
            }

            blst::G2Point out{};
            blst::G2Point::from_affine(out, affine);
            return out;
        }

        blst::GT deserialize_gt(const std::span<const std::uint8_t> bytes)
        {
            if (bytes.size() != (12 * 48)) {
                throw std::invalid_argument("BTX GT encoding has the wrong length");
            }

            blst::GT out{};
            blst::GT::deserialize(out, bytes.data());
            return out;
        }

        std::uint8_t read_u8(const std::span<const std::uint8_t> input, std::size_t& offset)
        {
            if (offset >= input.size()) {
                throw std::invalid_argument("BTX wire encoding is truncated");
            }
            return input[offset++];
        }

        std::uint32_t read_u32_be(const std::span<const std::uint8_t> input, std::size_t& offset)
        {
            if (input.size() - offset < 4) {
                throw std::invalid_argument("BTX wire encoding is truncated");
            }

            const std::uint32_t value =
                (static_cast<std::uint32_t>(input[offset]) << 24) |
                (static_cast<std::uint32_t>(input[offset + 1]) << 16) |
                (static_cast<std::uint32_t>(input[offset + 2]) << 8) |
                static_cast<std::uint32_t>(input[offset + 3]);
            offset += 4;
            return value;
        }

        std::span<const std::uint8_t> read_bytes(
            const std::span<const std::uint8_t> input,
            std::size_t& offset,
            const std::size_t count)
        {
            if (input.size() - offset < count) {
                throw std::invalid_argument("BTX wire encoding is truncated");
            }

            const std::span<const std::uint8_t> out = input.subspan(offset, count);
            offset += count;
            return out;
        }

        void expect_end(const std::span<const std::uint8_t> input, const std::size_t offset)
        {
            if (offset != input.size()) {
                throw std::invalid_argument("BTX wire encoding has trailing bytes");
            }
        }

        void expect_version(const std::span<const std::uint8_t> input, std::size_t& offset)
        {
            const std::uint8_t version = read_u8(input, offset);
            if (version != kBtxWireVersion) {
                throw std::invalid_argument("Unsupported BTX wire version");
            }
        }
    }

    Bytes SerializePublicParameters(const PublicParameters& public_parameters)
    {
        Bytes out;
        bte::building_blocks::append_u8(out, kBtxWireVersion);
        bte::building_blocks::append_u32_be(out, static_cast<std::uint32_t>(public_parameters.max_batch_size));
        bte::building_blocks::append_u32_be(out, static_cast<std::uint32_t>(public_parameters.num_servers));
        bte::building_blocks::append_u32_be(out, static_cast<std::uint32_t>(public_parameters.threshold));
        bte::building_blocks::append_array(out, serialize_gt(public_parameters.gt_tau_to_bmax_plus_one));
        bte::building_blocks::append_u32_be(out, static_cast<std::uint32_t>(public_parameters.g2_powers.size()));
        for (const blst::G2Point& point : public_parameters.g2_powers) {
            bte::building_blocks::append_array(out, compress_g2(point));
        }
        return out;
    }

    PublicParameters DeserializePublicParameters(const std::span<const std::uint8_t> encoded)
    {
        std::size_t offset = 0;
        expect_version(encoded, offset);

        PublicParameters public_parameters{};
        public_parameters.max_batch_size = read_u32_be(encoded, offset);
        public_parameters.num_servers = read_u32_be(encoded, offset);
        public_parameters.threshold = read_u32_be(encoded, offset);
        public_parameters.gt_tau_to_bmax_plus_one = deserialize_gt(read_bytes(encoded, offset, 12 * 48));

        const std::size_t num_powers = read_u32_be(encoded, offset);
        public_parameters.g2_powers.reserve(num_powers);
        for (std::size_t i = 0; i < num_powers; ++i) {
            public_parameters.g2_powers.push_back(deserialize_g2(read_bytes(encoded, offset, 96)));
        }

        expect_end(encoded, offset);
        return public_parameters;
    }

    Bytes SerializeCiphertext(const Ciphertext& ciphertext)
    {
        Bytes out;
        bte::building_blocks::append_u8(out, kBtxWireVersion);
        bte::building_blocks::append_array(out, compress_g1(ciphertext.c1));
        return out;
    }

    Ciphertext DeserializeCiphertext(const std::span<const std::uint8_t> encoded)
    {
        std::size_t offset = 0;
        expect_version(encoded, offset);

        Ciphertext ciphertext{};
        ciphertext.c1 = deserialize_g1(read_bytes(encoded, offset, 48));
        expect_end(encoded, offset);
        return ciphertext;
    }

    Bytes SerializeThresholdRequest(const ThresholdRequest& request)
    {
        Bytes out;
        bte::building_blocks::append_u8(out, kBtxWireVersion);
        bte::building_blocks::append_u32_be(out, static_cast<std::uint32_t>(request.ordered_c1.size()));
        for (const blst::G1Point& point : request.ordered_c1) {
            bte::building_blocks::append_array(out, compress_g1(point));
        }
        return out;
    }

    ThresholdRequest DeserializeThresholdRequest(const std::span<const std::uint8_t> encoded)
    {
        std::size_t offset = 0;
        expect_version(encoded, offset);

        ThresholdRequest request{};
        const std::size_t batch_size = read_u32_be(encoded, offset);
        request.ordered_c1.reserve(batch_size);
        for (std::size_t i = 0; i < batch_size; ++i) {
            request.ordered_c1.push_back(deserialize_g1(read_bytes(encoded, offset, 48)));
        }

        expect_end(encoded, offset);
        return request;
    }

    Bytes SerializePartialDecryptionShare(const PartialDecryptionShare& share)
    {
        Bytes out;
        bte::building_blocks::append_u8(out, kBtxWireVersion);
        bte::building_blocks::append_u32_be(out, share.index);
        bte::building_blocks::append_array(out, compress_g1(share.value));
        return out;
    }

    PartialDecryptionShare DeserializePartialDecryptionShare(const std::span<const std::uint8_t> encoded)
    {
        std::size_t offset = 0;
        expect_version(encoded, offset);

        PartialDecryptionShare share{};
        share.index = read_u32_be(encoded, offset);
        share.value = deserialize_g1(read_bytes(encoded, offset, 48));
        expect_end(encoded, offset);
        return share;
    }
}
