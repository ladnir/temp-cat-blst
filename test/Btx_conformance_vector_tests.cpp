#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "BTE/btx/admission.h"
#include "BTE/btx/core.h"

namespace
{
    enum class FieldKind
    {
        U64,
        Bytes,
        BytesList,
    };

    struct Field
    {
        FieldKind kind;
        std::string name;
        std::string encoded;
    };

    struct Case
    {
        std::string id;
        std::string operation;
        std::string expectation;
        std::vector<Field> fields;
    };

    bool valid_identifier(const std::string_view value)
    {
        if (value.empty()) {
            return false;
        }
        for (const char c : value) {
            if (!((c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_')) {
                return false;
            }
        }
        return true;
    }

    bool valid_bytes(const std::string_view value)
    {
        if (value == "-") {
            return true;
        }
        if (value.empty() || (value.size() & 1U) != 0U) {
            return false;
        }
        for (const char c : value) {
            if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'))) {
                return false;
            }
        }
        return true;
    }

    bool valid_bytes_list(const std::string_view value)
    {
        if (value == "-") {
            return true;
        }
        if (value.empty()) {
            return false;
        }
        std::size_t start = 0;
        for (;;) {
            const std::size_t comma = value.find(',', start);
            const std::string_view item = value.substr(
                start,
                comma == std::string_view::npos ? std::string_view::npos : comma - start);
            if (item == "-" || !valid_bytes(item)) {
                return false;
            }
            if (comma == std::string_view::npos) {
                return true;
            }
            start = comma + 1;
        }
    }

    bool valid_u64(const std::string_view value)
    {
        if (value.empty() || (value.size() > 1 && value.front() == '0')) {
            return false;
        }
        std::uint64_t parsed = 0;
        const char* const first = value.data();
        const char* const last = first + value.size();
        const auto result = std::from_chars(first, last, parsed);
        return result.ec == std::errc{} && result.ptr == last;
    }

    std::vector<std::string> split_words(const std::string& line)
    {
        std::istringstream stream(line);
        std::vector<std::string> words;
        for (std::string word; stream >> word;) {
            words.push_back(std::move(word));
        }
        return words;
    }

    std::vector<Case> parse_vectors(std::istream& input)
    {
        bool saw_header = false;
        bool in_case = false;
        Case current{};
        std::unordered_set<std::string> field_names;
        std::vector<Case> cases;

        std::string line;
        std::size_t line_number = 0;
        while (std::getline(input, line)) {
            ++line_number;
            const std::vector<std::string> words = split_words(line);
            if (words.empty() || words[0].starts_with('#')) {
                continue;
            }
            if (!saw_header) {
                if (words.size() != 1 || words[0] != "BTX-CONFORMANCE-V1") {
                    throw std::runtime_error("missing version header at line " + std::to_string(line_number));
                }
                saw_header = true;
                continue;
            }

            if (words.size() == 2 && words[0] == "case") {
                if (in_case || !valid_identifier(words[1])) {
                    throw std::runtime_error("invalid case at line " + std::to_string(line_number));
                }
                in_case = true;
                current = Case{};
                current.id = words[1];
                field_names.clear();
                continue;
            }
            if (words.size() == 2 && words[0] == "operation") {
                if (!in_case || !current.operation.empty() || !valid_identifier(words[1])) {
                    throw std::runtime_error("invalid operation at line " + std::to_string(line_number));
                }
                current.operation = words[1];
                continue;
            }
            if (words.size() == 2 && words[0] == "expect") {
                const bool valid_error = words[1].starts_with("error:") &&
                    valid_identifier(std::string_view(words[1]).substr(6));
                if (!in_case || !current.expectation.empty() || (words[1] != "ok" && !valid_error)) {
                    throw std::runtime_error("invalid expectation at line " + std::to_string(line_number));
                }
                current.expectation = words[1];
                continue;
            }
            if (words.size() == 3 &&
                (words[0] == "u64" || words[0] == "bytes" || words[0] == "bytes[]")) {
                if (!in_case || !valid_identifier(words[1]) || !field_names.insert(words[1]).second) {
                    throw std::runtime_error("invalid field at line " + std::to_string(line_number));
                }
                const FieldKind kind = words[0] == "u64"
                    ? FieldKind::U64
                    : (words[0] == "bytes" ? FieldKind::Bytes : FieldKind::BytesList);
                const bool valid = kind == FieldKind::U64
                    ? valid_u64(words[2])
                    : (kind == FieldKind::Bytes ? valid_bytes(words[2]) : valid_bytes_list(words[2]));
                if (!valid) {
                    throw std::runtime_error("invalid field value at line " + std::to_string(line_number));
                }
                current.fields.push_back(Field{ kind, words[1], words[2] });
                continue;
            }
            if (words.size() == 1 && words[0] == "end") {
                if (!in_case || current.operation.empty() || current.expectation.empty()) {
                    throw std::runtime_error("incomplete case at line " + std::to_string(line_number));
                }
                cases.push_back(std::move(current));
                current = Case{};
                field_names.clear();
                in_case = false;
                continue;
            }
            throw std::runtime_error("invalid syntax at line " + std::to_string(line_number));
        }

        if (!saw_header) {
            throw std::runtime_error("missing version header");
        }
        if (in_case) {
            throw std::runtime_error("unterminated case");
        }
        return cases;
    }

    const Field& find_field(const Case& test_case, const std::string_view name)
    {
        for (const Field& field : test_case.fields) {
            if (field.name == name) {
                return field;
            }
        }
        throw std::runtime_error("missing field " + std::string(name));
    }

    std::uint64_t u64_field(const Case& test_case, const std::string_view name)
    {
        const Field& field = find_field(test_case, name);
        if (field.kind != FieldKind::U64) {
            throw std::runtime_error("field is not u64");
        }
        std::uint64_t value = 0;
        const auto result = std::from_chars(
            field.encoded.data(),
            field.encoded.data() + field.encoded.size(),
            value);
        if (result.ec != std::errc{} || result.ptr != field.encoded.data() + field.encoded.size()) {
            throw std::runtime_error("invalid u64 field");
        }
        return value;
    }

    std::vector<std::uint8_t> bytes_field(const Case& test_case, const std::string_view name)
    {
        const Field& field = find_field(test_case, name);
        if (field.kind != FieldKind::Bytes) {
            throw std::runtime_error("field is not bytes");
        }
        if (field.encoded == "-") {
            return {};
        }

        const auto nibble = [](const char c) -> std::uint8_t {
            return c <= '9'
                ? static_cast<std::uint8_t>(c - '0')
                : static_cast<std::uint8_t>(c - 'a' + 10);
        };
        std::vector<std::uint8_t> bytes(field.encoded.size() / 2);
        for (std::size_t i = 0; i < bytes.size(); ++i) {
            bytes[i] = static_cast<std::uint8_t>(
                (nibble(field.encoded[2 * i]) << 4) | nibble(field.encoded[(2 * i) + 1]));
        }
        return bytes;
    }

    std::vector<std::uint32_t> decode_indices(const std::vector<std::uint8_t>& bytes)
    {
        if ((bytes.size() % 4) != 0) {
            throw std::runtime_error("index field is not a sequence of u32 values");
        }
        std::vector<std::uint32_t> indices;
        indices.reserve(bytes.size() / 4);
        for (std::size_t i = 0; i < bytes.size(); i += 4) {
            indices.push_back(
                (static_cast<std::uint32_t>(bytes[i]) << 24) |
                (static_cast<std::uint32_t>(bytes[i + 1]) << 16) |
                (static_cast<std::uint32_t>(bytes[i + 2]) << 8) |
                static_cast<std::uint32_t>(bytes[i + 3]));
        }
        return indices;
    }

    std::string combine_error_id(const std::invalid_argument& error)
    {
        const std::string_view message(error.what());
        if (message == "Need at least one BTX share to combine") {
            return "empty_shares";
        }
        if (message == "Need at least threshold-many BTX shares to combine") {
            return "insufficient_shares";
        }
        if (message == "BTX share indices must be nonzero") {
            return "zero_share_index";
        }
        if (message == "Duplicate validator index") {
            return "duplicate_share_index";
        }
        if (message == "BTX share index is outside the DKG domain") {
            return "invalid_share_index";
        }
        throw std::runtime_error("unmapped combination error: " + std::string(message));
    }
}

TEST(BtxConformanceVectors, ParserRejectsNonCanonicalInput)
{
    std::istringstream duplicate(
        "BTX-CONFORMANCE-V1\n"
        "case x\n"
        "operation format_smoke\n"
        "expect ok\n"
        "bytes a 00\n"
        "bytes a 01\n"
        "end\n");
    EXPECT_THROW(parse_vectors(duplicate), std::runtime_error);

    std::istringstream leading_zero(
        "BTX-CONFORMANCE-V1\n"
        "case x\n"
        "operation format_smoke\n"
        "expect ok\n"
        "u64 value 01\n"
        "end\n");
    EXPECT_THROW(parse_vectors(leading_zero), std::runtime_error);

    std::istringstream empty_list(
        "BTX-CONFORMANCE-V1\n"
        "case x\n"
        "operation format_smoke\n"
        "expect ok\n"
        "bytes[] values -\n"
        "end\n");
    EXPECT_EQ(parse_vectors(empty_list).size(), 1U);

    std::istringstream empty_element(
        "BTX-CONFORMANCE-V1\n"
        "case x\n"
        "operation format_smoke\n"
        "expect ok\n"
        "bytes[] values 00,-\n"
        "end\n");
    EXPECT_THROW(parse_vectors(empty_element), std::runtime_error);
}

TEST(BtxConformanceVectors, ParsesAuthoritativeSharedFile)
{
    const char* const path = std::getenv("BTX_CONFORMANCE_VECTOR_FILE");
    if (path == nullptr || path[0] == '\0') {
        GTEST_SKIP() << "BTX_CONFORMANCE_VECTOR_FILE is not configured";
    }

    std::ifstream input(path);
    ASSERT_TRUE(input.is_open()) << path;
    const std::vector<Case> cases = parse_vectors(input);
    ASSERT_EQ(cases.size(), 20U);
    EXPECT_EQ(cases[0].id, "combine_valid");
}

TEST(BtxConformanceVectors, CombineValidationMatchesAuthoritativeSharedFile)
{
    const char* const path = std::getenv("BTX_CONFORMANCE_VECTOR_FILE");
    if (path == nullptr || path[0] == '\0') {
        GTEST_SKIP() << "BTX_CONFORMANCE_VECTOR_FILE is not configured";
    }

    std::ifstream input(path);
    ASSERT_TRUE(input.is_open()) << path;
    const std::vector<Case> cases = parse_vectors(input);

    for (const Case& test_case : cases) {
        if (test_case.operation != "combine_validate") {
            continue;
        }

        const std::size_t threshold = static_cast<std::size_t>(
            u64_field(test_case, "reconstruction_threshold"));
        const std::size_t total_shares = static_cast<std::size_t>(
            u64_field(test_case, "total_shares"));
        const std::size_t share_count = static_cast<std::size_t>(
            u64_field(test_case, "share_count"));
        const std::vector<std::uint32_t> indices =
            decode_indices(bytes_field(test_case, "indices"));

        std::string actual;
        if (share_count != indices.size()) {
            actual = "error:mismatched_share_count";
        } else {
            bte::btx::PublicParameters public_parameters{};
            public_parameters.threshold = threshold;
            public_parameters.num_servers = total_shares;

            std::vector<bte::btx::PartialDecryptionShare> shares;
            shares.reserve(share_count);
            for (const std::uint32_t index : indices) {
                shares.push_back({ index, blst::G1Point::generator() });
            }

            try {
                (void)bte::btx::CombinePartialDecryptionShares(
                    std::span<const bte::btx::PartialDecryptionShare>(
                        shares.data(),
                        shares.size()),
                    public_parameters);
                actual = "ok";
            } catch (const std::invalid_argument& error) {
                actual = "error:" + combine_error_id(error);
            }
        }

        EXPECT_EQ(actual, test_case.expectation) << test_case.id;
    }
}

TEST(BtxConformanceVectors, CiphertextAdmissionMatchesAuthoritativeSharedFile)
{
    const char* const path = std::getenv("BTX_CONFORMANCE_VECTOR_FILE");
    if (path == nullptr || path[0] == '\0') {
        GTEST_SKIP() << "BTX_CONFORMANCE_VECTOR_FILE is not configured";
    }

    std::ifstream input(path);
    ASSERT_TRUE(input.is_open()) << path;
    const std::vector<Case> cases = parse_vectors(input);

    for (const Case& test_case : cases) {
        if (test_case.operation != "verify_ciphertext") {
            continue;
        }

        const std::vector<std::uint8_t> ciphertext =
            bytes_field(test_case, "ciphertext");
        const std::vector<std::uint8_t> associated_data =
            bytes_field(test_case, "ad");
        bte::btx::ClientCiphertext decoded{};
        const bte::btx::CiphertextAdmissionError error =
            bte::btx::DecodeAndVerifyClientCiphertext(
                ciphertext,
                associated_data,
                decoded);
        const std::string actual =
            error == bte::btx::CiphertextAdmissionError::Ok
            ? "ok"
            : "error:" + std::string(bte::btx::CiphertextAdmissionErrorId(error));
        EXPECT_EQ(actual, test_case.expectation) << test_case.id;
    }
}
