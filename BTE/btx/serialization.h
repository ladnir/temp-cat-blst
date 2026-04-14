#pragma once

#include <span>

#include "BTE/building_blocks/serialization.h"
#include "BTE/btx/types.h"

namespace bte::btx
{
    building_blocks::Bytes SerializePublicParameters(const PublicParameters& public_parameters);
    PublicParameters DeserializePublicParameters(std::span<const std::uint8_t> encoded);

    building_blocks::Bytes SerializeCiphertext(const Ciphertext& ciphertext);
    Ciphertext DeserializeCiphertext(std::span<const std::uint8_t> encoded);

    building_blocks::Bytes SerializeThresholdRequest(const ThresholdRequest& request);
    ThresholdRequest DeserializeThresholdRequest(std::span<const std::uint8_t> encoded);

    building_blocks::Bytes SerializePartialDecryptionShare(const PartialDecryptionShare& share);
    PartialDecryptionShare DeserializePartialDecryptionShare(std::span<const std::uint8_t> encoded);
}
