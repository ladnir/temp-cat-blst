#include <gtest/gtest.h>
#include <random>
#include <vector>
#include <iostream>

// Force both AVX512 and PORTABLE to be defined for this test file if possible,
// or check if they are defined in the build system.
// The user request says "These will be enabled if portable and avx512 are enabled."
// So we check for the defines.

#include "blst/common/defines.h"

// We need to inspect the internals of intrin_compat.h or ensure we can access both sets of functions.
// coping the logic from intrin_compat.h, if BLST_AVX512 is defined, we get intrinsics.
// if BLST_AVX512_PORTABLE is defined, we get the portable functions.

#if defined(BLST_AVX512) && defined(BLST_AVX512_PORTABLE)

#include "blst/backend/avx512/intrin_compat.h"

namespace blst {
namespace test {

// Helper to print mismatch
void print_failure(const char* op_name, int lane, u64 expected, u64 actual) {
    std::cerr << "Mismatch in " << op_name << " at lane " << lane << "\n"
              << "Expected: 0x" << std::hex << expected << "\n"
              << "Actual:   0x" << actual << std::dec << std::endl;
}

// Random number generator
std::mt19937_64 rng(12345);

cc512i random_cc512i() {
    u64 data[8];
    for (int i = 0; i < 8; ++i) data[i] = rng();
    return _cc512_loadu_si512(data); // Using the generic load which should work if available, or we construct manually
}

// Check equality between native intrinsic result and portable result
void check_equal(const cc512i& native_res, const cc512i& portable_res, const char* op_name) {
    // We need to extract values.
    // cc512i might be a wrapper around __m512i or a struct with u64 lane[8].
    // Since BLST_AVX512 is defined, cc512i usually has __m512i data; member if BLST_AVX_NATIVE is 1.
    // But BLST_AVX_NATIVE is defined inside intrin_compat.h based on logic.
    // If BLST_AVX512 is defined and BLST_AVX512_PORTABLE is defined, intrin_compat.h sets BLST_AVX_NATIVE to 0?
    // Let's re-read the file content logic provided earlier.
    
    /*
    #if defined(BLST_AVX512) && !defined(BLST_AVX512_PORTABLE)
    #define BLST_AVX_NATIVE 1
    ...
    #else
    #define BLST_AVX_NATIVE 0
    ...
    #endif
    */
    
    // Ah, if BOTH are defined, BLST_AVX_NATIVE is 0. This means the standard functions like `add_epi64`
    // map to `add_epi64_portable`. Use of intrinsics must be direct via `_mm512_*`.
    // Wait, if BLST_AVX_NATIVE is 0, then cc512i is defined as struct { u64 lane[8]; }.
    // If cc512i is u64 lane[8], we cannot pass it to _mm512_add_epi64 directly because that expects __m512i.
    
    // This is a problem for testing. The struct definition changes based on macros.
    // If we want to test compat, we need `cc512i` to be convertible or comparable.
    
    // If BLST_AVX_NATIVE is 0, cc512i is just lanes.
    // To test against intrinsics, we have to manually invoke intrinsics using casting or reinterpreting.
    
    u64 p_lanes[8];
    storeu_si512_portable(p_lanes, portable_res);
    
    u64 n_lanes[8];
    // If cc512i is the portable struct, native_res holds the portable result if we called the standard API.
    // We want to test "portable version against intrinsics".
    // But `check_equal` receives two cc512i's.
    // One should be computed by portable, one by native.
    // If cc512i is portable struct, we cannot compute native result easily into it using _mm512 functions
    // unless we store _mm512 result into memory and load it into cc512i.
    
    storeu_si512_portable(n_lanes, native_res);

    for (int i = 0; i < 8; ++i) {
        if (p_lanes[i] != n_lanes[i]) {
            print_failure(op_name, i, n_lanes[i], p_lanes[i]);
            ASSERT_EQ(p_lanes[i], n_lanes[i]) << op_name << " failed at lane " << i;
        }
    }
}

// Convert cc512i (which is likely portable struct here) to __m512i
__m512i to_native(const cc512i& a) {
    u64 tmp[8];
    storeu_si512_portable(tmp, a);
    return _mm512_loadu_si512(tmp);
}

// Convert __m512i to cc512i (portable struct)
cc512i from_native(__m512i a) {
    u64 tmp[8];
    _mm512_storeu_si512(tmp, a);
    return loadu_si512_portable(tmp);
}


class Avx512CompatTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Seed rng
    }
};

TEST_F(Avx512CompatTest, SetZero) {
    cc512i portable = setzero_si512_portable();
    cc512i native = from_native(_mm512_setzero_si512());
    check_equal(native, portable, "setzero_si512");
}

TEST_F(Avx512CompatTest, Set1Epi64) {
    for (int i = 0; i < 20; ++i) {
        u64 val = rng();
        cc512i portable = set1_epi64_portable(val);
        cc512i native = from_native(_mm512_set1_epi64(val));
        check_equal(native, portable, "set1_epi64");
    }
}

TEST_F(Avx512CompatTest, SetEpi64) {
    for (int i = 0; i < 20; ++i) {
        u64 x[8];
        for(int j=0; j<8; ++j) x[j] = rng();
        
        cc512i portable = set_epi64_portable(x[7], x[6], x[5], x[4], x[3], x[2], x[1], x[0]);
        cc512i native = from_native(_mm512_set_epi64(x[7], x[6], x[5], x[4], x[3], x[2], x[1], x[0]));
        check_equal(native, portable, "set_epi64");
    }
}

TEST_F(Avx512CompatTest, AddEpi64) {
    for (int i = 0; i < 100; ++i) {
        cc512i a = random_cc512i();
        cc512i b = random_cc512i();

        cc512i portable = add_epi64_portable(a, b);
        cc512i native = from_native(_mm512_add_epi64(to_native(a), to_native(b)));

        check_equal(native, portable, "add_epi64");
    }
}

TEST_F(Avx512CompatTest, SubEpi64) {
    for (int i = 0; i < 100; ++i) {
        cc512i a = random_cc512i();
        cc512i b = random_cc512i();

        cc512i portable = sub_epi64_portable(a, b);
        cc512i native = from_native(_mm512_sub_epi64(to_native(a), to_native(b)));

        check_equal(native, portable, "sub_epi64");
    }
}

TEST_F(Avx512CompatTest, LogicOps) {
    for (int i = 0; i < 100; ++i) {
        cc512i a = random_cc512i();
        cc512i b = random_cc512i();

        check_equal(from_native(_mm512_and_si512(to_native(a), to_native(b))), 
                    and_si512_portable(a, b), "and_si512");

        check_equal(from_native(_mm512_or_si512(to_native(a), to_native(b))), 
                    or_si512_portable(a, b), "or_si512");

        check_equal(from_native(_mm512_xor_si512(to_native(a), to_native(b))), 
                    xor_si512_portable(a, b), "xor_si512");
    }
}

TEST_F(Avx512CompatTest, Shifts) {
    for (int i = 0; i < 20; ++i) {
        cc512i a = random_cc512i();
        int bits = rng() % 64;

        check_equal(from_native(_mm512_slli_epi64(to_native(a), bits)),
                    slli_epi64_portable(a, bits), "slli_epi64");
        
        check_equal(from_native(_mm512_srli_epi64(to_native(a), bits)),
                    srli_epi64_portable(a, bits), "srli_epi64");

        check_equal(from_native(_mm512_srai_epi64(to_native(a), bits)),
                    srai_epi64_portable(a, bits), "srai_epi64");
    }
}

TEST_F(Avx512CompatTest, VariableShifts) {
    for (int i = 0; i < 20; ++i) {
        cc512i a = random_cc512i();
        cc512i bits = random_cc512i(); // Shifts only look at low bits but safe to randomize all

        check_equal(from_native(_mm512_sllv_epi64(to_native(a), to_native(bits))),
                    sllv_epi64_portable(a, bits), "sllv_epi64");
    }
}

TEST_F(Avx512CompatTest, Madd52) {
    for (int i = 0; i < 100; ++i) {
        cc512i acc = random_cc512i();
        cc512i a = random_cc512i();
        cc512i b = random_cc512i();

        check_equal(from_native(_mm512_madd52lo_epu64(to_native(acc), to_native(a), to_native(b))),
                    madd52lo_epu64_portable(acc, a, b), "madd52lo_epu64");

        check_equal(from_native(_mm512_madd52hi_epu64(to_native(acc), to_native(a), to_native(b))),
                    madd52hi_epu64_portable(acc, a, b), "madd52hi_epu64");
    }
}

TEST_F(Avx512CompatTest, MaskedOps) {
    for (int i = 0; i < 100; ++i) {
        cc512i src = random_cc512i();
        cc512i a = random_cc512i();
        cc512i b = random_cc512i();
        u64 mask = rng() & 0xFF; // 8-bit mask

        check_equal(from_native(_mm512_mask_add_epi64(to_native(src), (__mmask8)mask, to_native(a), to_native(b))),
                    mask_add_epi64_portable(src, mask, a, b), "mask_add_epi64");

        check_equal(from_native(_mm512_mask_sub_epi64(to_native(src), (__mmask8)mask, to_native(a), to_native(b))),
                    mask_sub_epi64_portable(src, mask, a, b), "mask_sub_epi64");
        
        check_equal(from_native(_mm512_mask_blend_epi64((__mmask8)mask, to_native(a), to_native(b))),
                    mask_blend_epi64_portable(mask, a, b), "mask_blend_epi64");
    }
}

TEST_F(Avx512CompatTest, Permutes) {
    for (int i = 0; i < 50; ++i) {
        cc512i a = random_cc512i();
        cc512i idx = random_cc512i();

        check_equal(from_native(_mm512_permutexvar_epi64(to_native(idx), to_native(a))),
                    permutexvar_epi64_portable(idx, a), "permutexvar_epi64");
    }
}

TEST_F(Avx512CompatTest, TemplatePermutes) {
    cc512i a = random_cc512i();
    
    // Test immediate permute/shuffle
    // Since these are templates we can just test a fixed meaningful immediate
    const int imm = 0b01001110; 

    // shuffle_epi32
    check_equal(from_native(_mm512_shuffle_epi32(to_native(a), (_MM_PERM_ENUM)imm)),
                shuffle_epi32_portable<imm>(a), "shuffle_epi32");
    
    // maskz_shuffle_epi32
    u64 mask = rng() & 0xFFFF; // 16-bit mask for epi32 (16 elements)
    check_equal(from_native(_mm512_maskz_shuffle_epi32((__mmask16)mask, to_native(a), (_MM_PERM_ENUM)imm)),
                maskz_shuffle_epi32_portable<imm>(mask, a), "maskz_shuffle_epi32");

    // permutex_epi64
    check_equal(from_native(_mm512_permutex_epi64(to_native(a), imm)),
                permutex_epi64_portable<imm>(a), "permutex_epi64");

    // maskz_permutex_epi64
    mask = rng() & 0xFF; // 8-bit mask for epi64
    check_equal(from_native(_mm512_maskz_permutex_epi64((__mmask8)mask, to_native(a), imm)),
                maskz_permutex_epi64_portable<imm>(mask, a), "maskz_permutex_epi64");
}

TEST_F(Avx512CompatTest, Extraction) {
    cc512i a = random_cc512i();

    // Extract 256-bit lower (imm8=0)
    cc256i low_portable = extracti64x4_epi64_portable<0>(a);
    cc256i low_native_portable_wrap{}; // Holds the result from native, but stored in portable struct if native defined
    
    // We need to convert cc256i as well
    // Native returns __m256i. We need to store it to compare.
    __m256i low_native = _mm512_extracti64x4_epi64(to_native(a), 0);
    
    u64 p_lanes[4];
    u64 n_lanes[4];
    
    // Manual store for 256 extraction test
    auto store256 = [](u64* dst, const cc256i& src) {
        // cc256i in portable mode has lane[4]
        // We can cast to u64* or access .lane
        for(int i=0; i<4; ++i) dst[i] = src.lane[i];
    };

    store256(p_lanes, low_portable);
    _mm256_storeu_si256((__m256i*)n_lanes, low_native);
    
    for(int i=0; i<4; ++i) ASSERT_EQ(p_lanes[i], n_lanes[i]) << "extracti64x4_epi64<0> mismatch lane " << i;

    // Extract 256-bit upper (imm8=1)
    cc256i high_portable = extracti64x4_epi64_portable<1>(a);
    __m256i high_native = _mm512_extracti64x4_epi64(to_native(a), 1);
    
    store256(p_lanes, high_portable);
    _mm256_storeu_si256((__m256i*)n_lanes, high_native);

    for(int i=0; i<4; ++i) ASSERT_EQ(p_lanes[i], n_lanes[i]) << "extracti64x4_epi64<1> mismatch lane " << i;

    // Extract 64-bit from 256-bit
    // Check all 4 indices
    EXPECT_EQ(extract_epi64_portable<0>(low_portable), (u64)_mm256_extract_epi64(low_native, 0));
    EXPECT_EQ(extract_epi64_portable<1>(low_portable), (u64)_mm256_extract_epi64(low_native, 1));
    EXPECT_EQ(extract_epi64_portable<2>(low_portable), (u64)_mm256_extract_epi64(low_native, 2));
    EXPECT_EQ(extract_epi64_portable<3>(low_portable), (u64)_mm256_extract_epi64(low_native, 3));
}

TEST_F(Avx512CompatTest, LoadStore) {
    u64 mem[8];
    for(int i=0; i<8; ++i) mem[i] = rng();

    // Load
    cc512i portable_loaded = loadu_si512_portable(mem);
    cc512i native_loaded = from_native(_mm512_loadu_si512(mem));
    check_equal(native_loaded, portable_loaded, "loadu_si512");
    
    // Store
    u64 mem_portable[8] = {0};
    u64 mem_native[8] = {0};
    
    storeu_si512_portable(mem_portable, portable_loaded);
    _mm512_storeu_si512(mem_native, to_native(portable_loaded));
    
    for(int i=0; i<8; ++i) ASSERT_EQ(mem_portable[i], mem_native[i]) << "storeu_si512 mismatch lane " << i;
    
    // Aligned load/store (simulated by just using aligned memory)
    alignas(64) u64 amem[8];
    for(int i=0; i<8; ++i) amem[i] = rng();
    
    // load_si512
    portable_loaded = load_si512_portable(amem);
    native_loaded = from_native(_mm512_load_si512(amem));
    check_equal(native_loaded, portable_loaded, "load_si512");
    
    // store_si512
    alignas(64) u64 amem_portable[8];
    alignas(64) u64 amem_native[8];
    
    store_si512_portable(amem_portable, portable_loaded);
    _mm512_store_si512(amem_native, to_native(portable_loaded));
    
    for(int i=0; i<8; ++i) ASSERT_EQ(amem_portable[i], amem_native[i]) << "store_si512 mismatch lane " << i;
}

} // namespace test
} // namespace blst

#endif // defined(BLST_AVX512) && defined(BLST_AVX512_PORTABLE)
