//
//#include <gtest/gtest.h>
//#include "../bindings/blst.hpp"
//namespace
//{
//
//
//    void printHex(const void* data, int byteSize)
//    {
//        unsigned char* byteData = (unsigned char*)data;
//        for (int i = 0; i < byteSize; i++) {
//            printf("%02x", byteData[i]);
//        }
//        printf("\n");
//    }
//}
//
//
//TEST(E1PointTest, CppApiCompressUncompress) {
//    using namespace blst;
//    {
//        const blst_p1* gen = blst_p1_generator();
//
//        uint8_t comp[48];
//        blst_p1_affine a1, a2;
//        blst_p1_to_affine(&a1, gen);
//        blst_p1_affine_compress(comp, &a1);
//
//        //printf("comp[0]=0x%02x\n", comp[0]);
//        //for (int i = 0; i < 48; i++)
//        //{
//        //    printf("%02x", comp[i]); puts("");
//        //}
//
//        BLST_ERROR err = blst_p1_uncompress(&a2, comp);
//        //printf("uncompress -> %d\n", err);
//        assert(err == BLST_SUCCESS);
//        assert(blst_p1_affine_on_curve(&a2) == 1);
//        assert(blst_p1_affine_in_g1(&a2) == 1);
//
//        blst_p1 p1, p2;
//        blst_p1_from_affine(&p1, &a1);
//        blst_p1_from_affine(&p2, &a2);
//        assert(blst_p1_is_equal(&p1, &p2) == 1);
//        //puts("OK");
//    }
//    {
//        uint8_t comp[48];
//        blst_p1_compress(comp, &*blst_p1_generator());
//
//        blst_p1_affine a;
//        ASSERT_EQ(blst_p1_uncompress(&a, comp), BLST_SUCCESS);
//
//        blst_p1 p;
//        blst_p1_from_affine(&p, &a);
//    }
//
//    // Create a point using the C++ API
//    P1 p1_cpp = P1::generator();
//
//    // Verify that the generator point is valid
//    EXPECT_TRUE(p1_cpp.on_curve());
//    EXPECT_TRUE(p1_cpp.in_group());
//
//    // Compress the point using the C++ API
//    std::vector<uint8_t> compressed(100);
//    p1_cpp.to_affine().compress(compressed.data());
//    //p1_cpp.compress(compressed.data());
//
//    //printf("C++ API compressed point (first byte: 0x%02X)\n", compressed[0]);
//	//printHex(compressed.data(), 48);
//
//    // Uncompress using the C++ API
//    P1 p2_cpp;
//    try {
//        p2_cpp = P1(compressed.data());
//
//        // Check the new point is valid
//        EXPECT_TRUE(p2_cpp.on_curve());
//        EXPECT_TRUE(p2_cpp.in_group());
//
//        // Compare with original point
//        EXPECT_TRUE(p1_cpp.is_equal(p2_cpp));
//
//        //printf("C++ API compression/decompression successful!\n");
//    }
//    catch (const std::exception& e) {
//        printf("C++ API decompression failed: %s\n", e.what());
//        FAIL() << "C++ API decompression failed";
//    }
//}
