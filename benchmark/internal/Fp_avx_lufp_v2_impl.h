#pragma once

#include <array>
#include <cstddef>
#include <cstdio>
#include <cstdlib>

#include "blst/fields/Fp.h"
#include "blst/fields/Fp12.h"
#include "blst/fields/Fp12_avx.h"
#include "blst/fields/Fp_avx.h"
#include "blst/fields/Fp6_avx.h"

namespace blst::bench_lufp_v2
{
#if BLST_AVX_NATIVE
    using v2_lane_mask_t = __mmask8;
#else
    using v2_lane_mask_t = u8;
#endif

    enum class LUFpV2RenormClass
    {
        ExactAdd,
        Fp6MulCommon,
    };

    inline constexpr std::array<u64, NWORDS> kV2P48X1 = {
        0xFFFFFFFFAAABULL, 0xB153FFFFB9FEULL, 0xF6241EABFFFEULL, 0x6730D2A0F6B0ULL,
        0x4B84F38512BFULL, 0x434BACD76477ULL, 0xE69A4B1BA7B6ULL, 0x1A0111EA397FULL,
    };

    inline const std::array<u64, NWORDS> kV2P48X2 = {
        0xFFFFFFFF5556ULL, 0x62A7FFFF73FDULL, 0xEC483D57FFFDULL, 0xCE61A541ED61ULL,
        0x9709E70A257EULL, 0x869759AEC8EEULL, 0xCD3496374F6CULL, 0x340223D472FFULL,
    };

    inline const std::array<u64, NWORDS> kV2P48X4 = {
        0xFFFFFFFEAAACULL, 0xC54FFFFEE7FBULL, 0xD8907AAFFFFAULL, 0x9CC34A83DAC3ULL,
        0x2E13CE144AFDULL, 0x0D2EB35D91DDULL, 0x9A692C6E9ED9ULL, 0x680447A8E5FFULL,
    };

    inline const std::array<u64, NWORDS> kV2P48X8 = {
        0xFFFFFFFD5558ULL, 0x8A9FFFFDCFF7ULL, 0xB120F55FFFF5ULL, 0x39869507B587ULL,
        0x5C279C2895FBULL, 0x1A5D66BB23BAULL, 0x34D258DD3DB2ULL, 0xD0088F51CBFFULL,
    };

    inline const std::array<u64, NWORDS> kV2P48X16 = {
        0xFFFFFFFAAAB0ULL, 0x153FFFFB9FEFULL, 0x6241EABFFFEBULL, 0x730D2A0F6B0FULL,
        0xB84F38512BF6ULL, 0x34BACD764774ULL, 0x69A4B1BA7B64ULL, 0x1A0111EA397FEULL,
    };

    inline const std::array<u64, NWORDS> kV2Bias3PPlus3 = {
        0xFFFFFFFF0004ULL, 0x13FBFFFF2DFCULL, 0xE26C5C03FFFCULL, 0x359277E2E412ULL,
        0xE28EDA8F383EULL, 0xC9E306862D65ULL, 0xB3CEE152F722ULL, 0x4E0335BEAC7FULL,
    };

    inline const std::array<cc512i, NWORDS> kV2P48X1Vec = {
        VSET(0x0000FFFFFFFFAAABLL, 0x0000FFFFFFFFAAABLL, 0x0000FFFFFFFFAAABLL, 0x0000FFFFFFFFAAABLL,
             0x0000FFFFFFFFAAABLL, 0x0000FFFFFFFFAAABLL, 0x0000FFFFFFFFAAABLL, 0x0000FFFFFFFFAAABLL),
        VSET(0x0000B153FFFFB9FELL, 0x0000B153FFFFB9FELL, 0x0000B153FFFFB9FELL, 0x0000B153FFFFB9FELL,
             0x0000B153FFFFB9FELL, 0x0000B153FFFFB9FELL, 0x0000B153FFFFB9FELL, 0x0000B153FFFFB9FELL),
        VSET(0x0000F6241EABFFFELL, 0x0000F6241EABFFFELL, 0x0000F6241EABFFFELL, 0x0000F6241EABFFFELL,
             0x0000F6241EABFFFELL, 0x0000F6241EABFFFELL, 0x0000F6241EABFFFELL, 0x0000F6241EABFFFELL),
        VSET(0x00006730D2A0F6B0LL, 0x00006730D2A0F6B0LL, 0x00006730D2A0F6B0LL, 0x00006730D2A0F6B0LL,
             0x00006730D2A0F6B0LL, 0x00006730D2A0F6B0LL, 0x00006730D2A0F6B0LL, 0x00006730D2A0F6B0LL),
        VSET(0x00004B84F38512BFLL, 0x00004B84F38512BFLL, 0x00004B84F38512BFLL, 0x00004B84F38512BFLL,
             0x00004B84F38512BFLL, 0x00004B84F38512BFLL, 0x00004B84F38512BFLL, 0x00004B84F38512BFLL),
        VSET(0x0000434BACD76477LL, 0x0000434BACD76477LL, 0x0000434BACD76477LL, 0x0000434BACD76477LL,
             0x0000434BACD76477LL, 0x0000434BACD76477LL, 0x0000434BACD76477LL, 0x0000434BACD76477LL),
        VSET(0x0000E69A4B1BA7B6LL, 0x0000E69A4B1BA7B6LL, 0x0000E69A4B1BA7B6LL, 0x0000E69A4B1BA7B6LL,
             0x0000E69A4B1BA7B6LL, 0x0000E69A4B1BA7B6LL, 0x0000E69A4B1BA7B6LL, 0x0000E69A4B1BA7B6LL),
        VSET(0x00001A0111EA397FLL, 0x00001A0111EA397FLL, 0x00001A0111EA397FLL, 0x00001A0111EA397FLL,
             0x00001A0111EA397FLL, 0x00001A0111EA397FLL, 0x00001A0111EA397FLL, 0x00001A0111EA397FLL),
    };

    inline const std::array<cc512i, NWORDS> kV2P48X2Vec = {
        VSET(0x0000FFFFFFFF5556LL, 0x0000FFFFFFFF5556LL, 0x0000FFFFFFFF5556LL, 0x0000FFFFFFFF5556LL,
             0x0000FFFFFFFF5556LL, 0x0000FFFFFFFF5556LL, 0x0000FFFFFFFF5556LL, 0x0000FFFFFFFF5556LL),
        VSET(0x000062A7FFFF73FDLL, 0x000062A7FFFF73FDLL, 0x000062A7FFFF73FDLL, 0x000062A7FFFF73FDLL,
             0x000062A7FFFF73FDLL, 0x000062A7FFFF73FDLL, 0x000062A7FFFF73FDLL, 0x000062A7FFFF73FDLL),
        VSET(0x0000EC483D57FFFDLL, 0x0000EC483D57FFFDLL, 0x0000EC483D57FFFDLL, 0x0000EC483D57FFFDLL,
             0x0000EC483D57FFFDLL, 0x0000EC483D57FFFDLL, 0x0000EC483D57FFFDLL, 0x0000EC483D57FFFDLL),
        VSET(0x0000CE61A541ED61LL, 0x0000CE61A541ED61LL, 0x0000CE61A541ED61LL, 0x0000CE61A541ED61LL,
             0x0000CE61A541ED61LL, 0x0000CE61A541ED61LL, 0x0000CE61A541ED61LL, 0x0000CE61A541ED61LL),
        VSET(0x00009709E70A257ELL, 0x00009709E70A257ELL, 0x00009709E70A257ELL, 0x00009709E70A257ELL,
             0x00009709E70A257ELL, 0x00009709E70A257ELL, 0x00009709E70A257ELL, 0x00009709E70A257ELL),
        VSET(0x0000869759AEC8EELL, 0x0000869759AEC8EELL, 0x0000869759AEC8EELL, 0x0000869759AEC8EELL,
             0x0000869759AEC8EELL, 0x0000869759AEC8EELL, 0x0000869759AEC8EELL, 0x0000869759AEC8EELL),
        VSET(0x0000CD3496374F6CLL, 0x0000CD3496374F6CLL, 0x0000CD3496374F6CLL, 0x0000CD3496374F6CLL,
             0x0000CD3496374F6CLL, 0x0000CD3496374F6CLL, 0x0000CD3496374F6CLL, 0x0000CD3496374F6CLL),
        VSET(0x0000340223D472FFLL, 0x0000340223D472FFLL, 0x0000340223D472FFLL, 0x0000340223D472FFLL,
             0x0000340223D472FFLL, 0x0000340223D472FFLL, 0x0000340223D472FFLL, 0x0000340223D472FFLL),
    };

    inline const std::array<u64, NWORDS> kV2P48X3 = {
        0xFFFFFFFF0001ULL, 0x13FBFFFF2DFCULL, 0xE26C5C03FFFCULL, 0x359277E2E412ULL,
        0xE28EDA8F383EULL, 0xC9E306862D65ULL, 0xB3CEE152F722ULL, 0x4E0335BEAC7FULL,
    };

    inline const std::array<cc512i, NWORDS> kV2P48X3Vec = {
        VSET(0x0000FFFFFFFF0001LL, 0x0000FFFFFFFF0001LL, 0x0000FFFFFFFF0001LL, 0x0000FFFFFFFF0001LL,
             0x0000FFFFFFFF0001LL, 0x0000FFFFFFFF0001LL, 0x0000FFFFFFFF0001LL, 0x0000FFFFFFFF0001LL),
        VSET(0x000013FBFFFF2DFCLL, 0x000013FBFFFF2DFCLL, 0x000013FBFFFF2DFCLL, 0x000013FBFFFF2DFCLL,
             0x000013FBFFFF2DFCLL, 0x000013FBFFFF2DFCLL, 0x000013FBFFFF2DFCLL, 0x000013FBFFFF2DFCLL),
        VSET(0x0000E26C5C03FFFCLL, 0x0000E26C5C03FFFCLL, 0x0000E26C5C03FFFCLL, 0x0000E26C5C03FFFCLL,
             0x0000E26C5C03FFFCLL, 0x0000E26C5C03FFFCLL, 0x0000E26C5C03FFFCLL, 0x0000E26C5C03FFFCLL),
        VSET(0x0000359277E2E412LL, 0x0000359277E2E412LL, 0x0000359277E2E412LL, 0x0000359277E2E412LL,
             0x0000359277E2E412LL, 0x0000359277E2E412LL, 0x0000359277E2E412LL, 0x0000359277E2E412LL),
        VSET(0x0000E28EDA8F383ELL, 0x0000E28EDA8F383ELL, 0x0000E28EDA8F383ELL, 0x0000E28EDA8F383ELL,
             0x0000E28EDA8F383ELL, 0x0000E28EDA8F383ELL, 0x0000E28EDA8F383ELL, 0x0000E28EDA8F383ELL),
        VSET(0x0000C9E306862D65LL, 0x0000C9E306862D65LL, 0x0000C9E306862D65LL, 0x0000C9E306862D65LL,
             0x0000C9E306862D65LL, 0x0000C9E306862D65LL, 0x0000C9E306862D65LL, 0x0000C9E306862D65LL),
        VSET(0x0000B3CEE152F722LL, 0x0000B3CEE152F722LL, 0x0000B3CEE152F722LL, 0x0000B3CEE152F722LL,
             0x0000B3CEE152F722LL, 0x0000B3CEE152F722LL, 0x0000B3CEE152F722LL, 0x0000B3CEE152F722LL),
        VSET(0x00004E0335BEAC7FLL, 0x00004E0335BEAC7FLL, 0x00004E0335BEAC7FLL, 0x00004E0335BEAC7FLL,
             0x00004E0335BEAC7FLL, 0x00004E0335BEAC7FLL, 0x00004E0335BEAC7FLL, 0x00004E0335BEAC7FLL),
    };

    inline const std::array<cc512i, NWORDS> kV2P48X4Vec = {
        VSET(0x0000FFFFFFFEAAACLL, 0x0000FFFFFFFEAAACLL, 0x0000FFFFFFFEAAACLL, 0x0000FFFFFFFEAAACLL,
             0x0000FFFFFFFEAAACLL, 0x0000FFFFFFFEAAACLL, 0x0000FFFFFFFEAAACLL, 0x0000FFFFFFFEAAACLL),
        VSET(0x0000C54FFFFEE7FBLL, 0x0000C54FFFFEE7FBLL, 0x0000C54FFFFEE7FBLL, 0x0000C54FFFFEE7FBLL,
             0x0000C54FFFFEE7FBLL, 0x0000C54FFFFEE7FBLL, 0x0000C54FFFFEE7FBLL, 0x0000C54FFFFEE7FBLL),
        VSET(0x0000D8907AAFFFFALL, 0x0000D8907AAFFFFALL, 0x0000D8907AAFFFFALL, 0x0000D8907AAFFFFALL,
             0x0000D8907AAFFFFALL, 0x0000D8907AAFFFFALL, 0x0000D8907AAFFFFALL, 0x0000D8907AAFFFFALL),
        VSET(0x00009CC34A83DAC3LL, 0x00009CC34A83DAC3LL, 0x00009CC34A83DAC3LL, 0x00009CC34A83DAC3LL,
             0x00009CC34A83DAC3LL, 0x00009CC34A83DAC3LL, 0x00009CC34A83DAC3LL, 0x00009CC34A83DAC3LL),
        VSET(0x00002E13CE144AFDLL, 0x00002E13CE144AFDLL, 0x00002E13CE144AFDLL, 0x00002E13CE144AFDLL,
             0x00002E13CE144AFDLL, 0x00002E13CE144AFDLL, 0x00002E13CE144AFDLL, 0x00002E13CE144AFDLL),
        VSET(0x00000D2EB35D91DDLL, 0x00000D2EB35D91DDLL, 0x00000D2EB35D91DDLL, 0x00000D2EB35D91DDLL,
             0x00000D2EB35D91DDLL, 0x00000D2EB35D91DDLL, 0x00000D2EB35D91DDLL, 0x00000D2EB35D91DDLL),
        VSET(0x00009A692C6E9ED9LL, 0x00009A692C6E9ED9LL, 0x00009A692C6E9ED9LL, 0x00009A692C6E9ED9LL,
             0x00009A692C6E9ED9LL, 0x00009A692C6E9ED9LL, 0x00009A692C6E9ED9LL, 0x00009A692C6E9ED9LL),
        VSET(0x0000680447A8E5FFLL, 0x0000680447A8E5FFLL, 0x0000680447A8E5FFLL, 0x0000680447A8E5FFLL,
             0x0000680447A8E5FFLL, 0x0000680447A8E5FFLL, 0x0000680447A8E5FFLL, 0x0000680447A8E5FFLL),
    };

    inline const std::array<cc512i, NWORDS> kV2P48X8Vec = {
        VSET(0x0000FFFFFFFD5558LL, 0x0000FFFFFFFD5558LL, 0x0000FFFFFFFD5558LL, 0x0000FFFFFFFD5558LL,
             0x0000FFFFFFFD5558LL, 0x0000FFFFFFFD5558LL, 0x0000FFFFFFFD5558LL, 0x0000FFFFFFFD5558LL),
        VSET(0x00008A9FFFFDCFF7LL, 0x00008A9FFFFDCFF7LL, 0x00008A9FFFFDCFF7LL, 0x00008A9FFFFDCFF7LL,
             0x00008A9FFFFDCFF7LL, 0x00008A9FFFFDCFF7LL, 0x00008A9FFFFDCFF7LL, 0x00008A9FFFFDCFF7LL),
        VSET(0x0000B120F55FFFF5LL, 0x0000B120F55FFFF5LL, 0x0000B120F55FFFF5LL, 0x0000B120F55FFFF5LL,
             0x0000B120F55FFFF5LL, 0x0000B120F55FFFF5LL, 0x0000B120F55FFFF5LL, 0x0000B120F55FFFF5LL),
        VSET(0x000039869507B587LL, 0x000039869507B587LL, 0x000039869507B587LL, 0x000039869507B587LL,
             0x000039869507B587LL, 0x000039869507B587LL, 0x000039869507B587LL, 0x000039869507B587LL),
        VSET(0x00005C279C2895FBLL, 0x00005C279C2895FBLL, 0x00005C279C2895FBLL, 0x00005C279C2895FBLL,
             0x00005C279C2895FBLL, 0x00005C279C2895FBLL, 0x00005C279C2895FBLL, 0x00005C279C2895FBLL),
        VSET(0x00001A5D66BB23BALL, 0x00001A5D66BB23BALL, 0x00001A5D66BB23BALL, 0x00001A5D66BB23BALL,
             0x00001A5D66BB23BALL, 0x00001A5D66BB23BALL, 0x00001A5D66BB23BALL, 0x00001A5D66BB23BALL),
        VSET(0x000034D258DD3DB2LL, 0x000034D258DD3DB2LL, 0x000034D258DD3DB2LL, 0x000034D258DD3DB2LL,
             0x000034D258DD3DB2LL, 0x000034D258DD3DB2LL, 0x000034D258DD3DB2LL, 0x000034D258DD3DB2LL),
        VSET(0x0000D0088F51CBFFLL, 0x0000D0088F51CBFFLL, 0x0000D0088F51CBFFLL, 0x0000D0088F51CBFFLL,
             0x0000D0088F51CBFFLL, 0x0000D0088F51CBFFLL, 0x0000D0088F51CBFFLL, 0x0000D0088F51CBFFLL),
    };

    inline const std::array<cc512i, NWORDS> kV2P48X16Vec = {
        VSET(0x0000FFFFFFFAAAB0LL, 0x0000FFFFFFFAAAB0LL, 0x0000FFFFFFFAAAB0LL, 0x0000FFFFFFFAAAB0LL,
             0x0000FFFFFFFAAAB0LL, 0x0000FFFFFFFAAAB0LL, 0x0000FFFFFFFAAAB0LL, 0x0000FFFFFFFAAAB0LL),
        VSET(0x0000153FFFFB9FEFLL, 0x0000153FFFFB9FEFLL, 0x0000153FFFFB9FEFLL, 0x0000153FFFFB9FEFLL,
             0x0000153FFFFB9FEFLL, 0x0000153FFFFB9FEFLL, 0x0000153FFFFB9FEFLL, 0x0000153FFFFB9FEFLL),
        VSET(0x00006241EABFFFEBLL, 0x00006241EABFFFEBLL, 0x00006241EABFFFEBLL, 0x00006241EABFFFEBLL,
             0x00006241EABFFFEBLL, 0x00006241EABFFFEBLL, 0x00006241EABFFFEBLL, 0x00006241EABFFFEBLL),
        VSET(0x0000730D2A0F6B0FLL, 0x0000730D2A0F6B0FLL, 0x0000730D2A0F6B0FLL, 0x0000730D2A0F6B0FLL,
             0x0000730D2A0F6B0FLL, 0x0000730D2A0F6B0FLL, 0x0000730D2A0F6B0FLL, 0x0000730D2A0F6B0FLL),
        VSET(0x0000B84F38512BF6LL, 0x0000B84F38512BF6LL, 0x0000B84F38512BF6LL, 0x0000B84F38512BF6LL,
             0x0000B84F38512BF6LL, 0x0000B84F38512BF6LL, 0x0000B84F38512BF6LL, 0x0000B84F38512BF6LL),
        VSET(0x000034BACD764774LL, 0x000034BACD764774LL, 0x000034BACD764774LL, 0x000034BACD764774LL,
             0x000034BACD764774LL, 0x000034BACD764774LL, 0x000034BACD764774LL, 0x000034BACD764774LL),
        VSET(0x000069A4B1BA7B64LL, 0x000069A4B1BA7B64LL, 0x000069A4B1BA7B64LL, 0x000069A4B1BA7B64LL,
             0x000069A4B1BA7B64LL, 0x000069A4B1BA7B64LL, 0x000069A4B1BA7B64LL, 0x000069A4B1BA7B64LL),
        VSET(0x0001A0111EA397FELL, 0x0001A0111EA397FELL, 0x0001A0111EA397FELL, 0x0001A0111EA397FELL,
             0x0001A0111EA397FELL, 0x0001A0111EA397FELL, 0x0001A0111EA397FELL, 0x0001A0111EA397FELL),
    };

    inline const std::array<cc512i, NWORDS> kV2Bias3PPlus3Vec = {
        VSET(0x0000FFFFFFFF0004LL, 0x0000FFFFFFFF0004LL, 0x0000FFFFFFFF0004LL, 0x0000FFFFFFFF0004LL,
             0x0000FFFFFFFF0004LL, 0x0000FFFFFFFF0004LL, 0x0000FFFFFFFF0004LL, 0x0000FFFFFFFF0004LL),
        VSET(0x000013FBFFFF2DFCLL, 0x000013FBFFFF2DFCLL, 0x000013FBFFFF2DFCLL, 0x000013FBFFFF2DFCLL,
             0x000013FBFFFF2DFCLL, 0x000013FBFFFF2DFCLL, 0x000013FBFFFF2DFCLL, 0x000013FBFFFF2DFCLL),
        VSET(0x0000E26C5C03FFFCLL, 0x0000E26C5C03FFFCLL, 0x0000E26C5C03FFFCLL, 0x0000E26C5C03FFFCLL,
             0x0000E26C5C03FFFCLL, 0x0000E26C5C03FFFCLL, 0x0000E26C5C03FFFCLL, 0x0000E26C5C03FFFCLL),
        VSET(0x0000359277E2E412LL, 0x0000359277E2E412LL, 0x0000359277E2E412LL, 0x0000359277E2E412LL,
             0x0000359277E2E412LL, 0x0000359277E2E412LL, 0x0000359277E2E412LL, 0x0000359277E2E412LL),
        VSET(0x0000E28EDA8F383ELL, 0x0000E28EDA8F383ELL, 0x0000E28EDA8F383ELL, 0x0000E28EDA8F383ELL,
             0x0000E28EDA8F383ELL, 0x0000E28EDA8F383ELL, 0x0000E28EDA8F383ELL, 0x0000E28EDA8F383ELL),
        VSET(0x0000C9E306862D65LL, 0x0000C9E306862D65LL, 0x0000C9E306862D65LL, 0x0000C9E306862D65LL,
             0x0000C9E306862D65LL, 0x0000C9E306862D65LL, 0x0000C9E306862D65LL, 0x0000C9E306862D65LL),
        VSET(0x0000B3CEE152F722LL, 0x0000B3CEE152F722LL, 0x0000B3CEE152F722LL, 0x0000B3CEE152F722LL,
             0x0000B3CEE152F722LL, 0x0000B3CEE152F722LL, 0x0000B3CEE152F722LL, 0x0000B3CEE152F722LL),
        VSET(0x00004E0335BEAC7FLL, 0x00004E0335BEAC7FLL, 0x00004E0335BEAC7FLL, 0x00004E0335BEAC7FLL,
             0x00004E0335BEAC7FLL, 0x00004E0335BEAC7FLL, 0x00004E0335BEAC7FLL, 0x00004E0335BEAC7FLL),
    };

    template <u64 Scale>
    BLST_FORCEINLINE constexpr auto v2_make_scaled_words() -> std::array<u64, NWORDS>
    {
        std::array<u64, NWORDS> out{};
        u64 carry = 0;

        for (std::size_t limb = 0; limb + 1 < NWORDS; ++limb) {
            const u64 prod = kV2P48X1[limb] * Scale + carry;
            out[limb] = prod & BMASK;
            carry = prod >> BRADIX;
        }

        out[NWORDS - 1] = kV2P48X1[NWORDS - 1] * Scale + carry;
        return out;
    }

    BLST_FORCEINLINE auto v2_make_scaled_vecs(const std::array<u64, NWORDS>& words) -> std::array<cc512i, NWORDS>
    {
        std::array<cc512i, NWORDS> out{};
        for (std::size_t limb = 0; limb < NWORDS; ++limb) {
            out[limb] = VSET1(static_cast<i64>(words[limb]));
        }
        return out;
    }

    inline constexpr auto kV2P48X5 = v2_make_scaled_words<5>();
    inline constexpr auto kV2P48X6 = v2_make_scaled_words<6>();
    inline constexpr auto kV2P48X7 = v2_make_scaled_words<7>();
    inline constexpr auto kV2P48X9 = v2_make_scaled_words<9>();
    inline constexpr auto kV2P48X10 = v2_make_scaled_words<10>();
    inline constexpr auto kV2P48X11 = v2_make_scaled_words<11>();
    inline constexpr auto kV2P48X12 = v2_make_scaled_words<12>();
    inline constexpr auto kV2P48X13 = v2_make_scaled_words<13>();
    inline constexpr auto kV2P48X14 = v2_make_scaled_words<14>();
    inline constexpr auto kV2P48X15 = v2_make_scaled_words<15>();

    inline const auto kV2P48X5Vec = v2_make_scaled_vecs(kV2P48X5);
    inline const auto kV2P48X6Vec = v2_make_scaled_vecs(kV2P48X6);
    inline const auto kV2P48X7Vec = v2_make_scaled_vecs(kV2P48X7);
    inline const auto kV2P48X9Vec = v2_make_scaled_vecs(kV2P48X9);
    inline const auto kV2P48X10Vec = v2_make_scaled_vecs(kV2P48X10);
    inline const auto kV2P48X11Vec = v2_make_scaled_vecs(kV2P48X11);
    inline const auto kV2P48X12Vec = v2_make_scaled_vecs(kV2P48X12);
    inline const auto kV2P48X13Vec = v2_make_scaled_vecs(kV2P48X13);
    inline const auto kV2P48X14Vec = v2_make_scaled_vecs(kV2P48X14);
    inline const auto kV2P48X15Vec = v2_make_scaled_vecs(kV2P48X15);

    // Returns the scalar radix-48 representation of Scale * p for the supported small scales.
    template <u64 Scale>
    BLST_FORCEINLINE constexpr auto v2_scaled_words() -> const std::array<u64, NWORDS>&
    {
        if constexpr (Scale == 1) {
            return kV2P48X1;
        }
        else if constexpr (Scale == 2) {
            return kV2P48X2;
        }
        else if constexpr (Scale == 3) {
            return kV2P48X3;
        }
        else if constexpr (Scale == 4) {
            return kV2P48X4;
        }
        else if constexpr (Scale == 5) {
            return kV2P48X5;
        }
        else if constexpr (Scale == 6) {
            return kV2P48X6;
        }
        else if constexpr (Scale == 7) {
            return kV2P48X7;
        }
        else if constexpr (Scale == 8) {
            return kV2P48X8;
        }
        else if constexpr (Scale == 9) {
            return kV2P48X9;
        }
        else if constexpr (Scale == 10) {
            return kV2P48X10;
        }
        else if constexpr (Scale == 11) {
            return kV2P48X11;
        }
        else if constexpr (Scale == 12) {
            return kV2P48X12;
        }
        else if constexpr (Scale == 13) {
            return kV2P48X13;
        }
        else if constexpr (Scale == 14) {
            return kV2P48X14;
        }
        else if constexpr (Scale == 15) {
            return kV2P48X15;
        }
        else if constexpr (Scale == 16) {
            return kV2P48X16;
        }
        else {
            static_assert(
                (Scale >= 1 && Scale <= 16),
                "unsupported v2 sign-step scale");
            return kV2P48X1;
        }
    }

    // Returns the vectorized radix-48 representation of Scale * p for the supported small scales.
    template <u64 Scale>
    BLST_FORCEINLINE auto v2_scaled_vecs() -> const std::array<cc512i, NWORDS>&
    {
        if constexpr (Scale == 1) {
            return kV2P48X1Vec;
        }
        else if constexpr (Scale == 2) {
            return kV2P48X2Vec;
        }
        else if constexpr (Scale == 3) {
            return kV2P48X3Vec;
        }
        else if constexpr (Scale == 4) {
            return kV2P48X4Vec;
        }
        else if constexpr (Scale == 5) {
            return kV2P48X5Vec;
        }
        else if constexpr (Scale == 6) {
            return kV2P48X6Vec;
        }
        else if constexpr (Scale == 7) {
            return kV2P48X7Vec;
        }
        else if constexpr (Scale == 8) {
            return kV2P48X8Vec;
        }
        else if constexpr (Scale == 9) {
            return kV2P48X9Vec;
        }
        else if constexpr (Scale == 10) {
            return kV2P48X10Vec;
        }
        else if constexpr (Scale == 11) {
            return kV2P48X11Vec;
        }
        else if constexpr (Scale == 12) {
            return kV2P48X12Vec;
        }
        else if constexpr (Scale == 13) {
            return kV2P48X13Vec;
        }
        else if constexpr (Scale == 14) {
            return kV2P48X14Vec;
        }
        else if constexpr (Scale == 15) {
            return kV2P48X15Vec;
        }
        else if constexpr (Scale == 16) {
            return kV2P48X16Vec;
        }
        else {
            static_assert(
                (Scale >= 1 && Scale <= 16),
                "unsupported v2 sign-step scale");
            return kV2P48X1Vec;
        }
    }

    // Lane-wise signed a > b compare that matches the current AVX/native portability layer.
    BLST_FORCEINLINE auto v2_mask_gt_epi64(const cc512i a, const cc512i b) -> v2_lane_mask_t
    {
#if BLST_AVX_NATIVE
        return _mm512_cmpgt_epi64_mask(to_native_512(a), to_native_512(b));
#else
        std::array<u64, 8> lhs{};
        std::array<u64, 8> rhs{};
        _cc512_storeu_si512(lhs.data(), a);
        _cc512_storeu_si512(rhs.data(), b);
        u8 mask = 0;
        for (u8 lane = 0; lane < 8; ++lane) {
            if (static_cast<i64>(lhs[lane]) > static_cast<i64>(rhs[lane])) {
                mask |= static_cast<u8>(1u << lane);
            }
        }
        return mask;
#endif
    }

    // Lane-wise a == b compare that returns the native lane-mask type.
    BLST_FORCEINLINE auto v2_mask_eq_epi64(const cc512i a, const cc512i b) -> v2_lane_mask_t
    {
#if BLST_AVX_NATIVE
        return _mm512_cmpeq_epi64_mask(to_native_512(a), to_native_512(b));
#else
        std::array<u64, 8> lhs{};
        std::array<u64, 8> rhs{};
        _cc512_storeu_si512(lhs.data(), a);
        _cc512_storeu_si512(rhs.data(), b);
        u8 mask = 0;
        for (u8 lane = 0; lane < 8; ++lane) {
            if (lhs[lane] == rhs[lane]) {
                mask |= static_cast<u8>(1u << lane);
            }
        }
        return mask;
#endif
    }

    // Lane-wise a == b compare collapsed to the benchmark u8 mask format.
    BLST_FORCEINLINE auto v2_mask_eq_epi64_u8(const cc512i a, const cc512i b) -> u8
    {
#if BLST_AVX_NATIVE
        return static_cast<u8>(_mm512_cmpeq_epi64_mask(to_native_512(a), to_native_512(b)));
#else
        std::array<u64, 8> lhs{};
        std::array<u64, 8> rhs{};
        _cc512_storeu_si512(lhs.data(), a);
        _cc512_storeu_si512(rhs.data(), b);
        u8 mask = 0;
        for (u8 lane = 0; lane < 8; ++lane) {
            if (lhs[lane] == rhs[lane]) {
                mask |= static_cast<u8>(1u << lane);
            }
        }
        return mask;
#endif
    }

    // Lane-wise signed a > b compare collapsed to the benchmark u8 mask format.
    BLST_FORCEINLINE auto v2_mask_gt_epi64_u8(const cc512i a, const cc512i b) -> u8
    {
#if BLST_AVX_NATIVE
        return static_cast<u8>(_mm512_cmpgt_epi64_mask(to_native_512(a), to_native_512(b)));
#else
        std::array<u64, 8> lhs{};
        std::array<u64, 8> rhs{};
        _cc512_storeu_si512(lhs.data(), a);
        _cc512_storeu_si512(rhs.data(), b);
        u8 mask = 0;
        for (u8 lane = 0; lane < 8; ++lane) {
            if (static_cast<i64>(lhs[lane]) > static_cast<i64>(rhs[lane])) {
                mask |= static_cast<u8>(1u << lane);
            }
        }
        return mask;
#endif
    }

    // Extracts the sign bit of each lane as a lane mask.
    BLST_FORCEINLINE auto v2_signbit_mask(const cc512i a) -> v2_lane_mask_t
    {
#if BLST_AVX_NATIVE
        return _mm512_movepi64_mask(to_native_512(a));
#else
        std::array<u64, 8> lanes{};
        _cc512_storeu_si512(lanes.data(), a);
        u8 mask = 0;
        for (u8 lane = 0; lane < 8; ++lane) {
            mask |= static_cast<u8>(((lanes[lane] >> 63) & 1u) << lane);
        }
        return mask;
#endif
    }
    //
    //BLST_FORCEINLINE auto lufp_v2_alias(UFpPack8x1& src) -> UFpPack8x1&
    //{
    //    return reinterpret_cast<UFpPack8x1&>(src);
    //}
    //
    //BLST_FORCEINLINE auto lufp_v2_alias(const UFpPack8x1& src) -> const UFpPack8x1&
    //{
    //    return reinterpret_cast<const UFpPack8x1&>(src);
    //}

    // Limb-wise widened add without normalization.
    BLST_FORCEINLINE void lufp_v2_add(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
    {
        for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
            r[limb] = VADD(a[limb], b[limb]);
        }
    }

    // Limb-wise widened subtract without normalization.
    BLST_FORCEINLINE void lufp_v2_sub(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
    {
        for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
            r[limb] = VSUB(a[limb], b[limb]);
        }
    }

    // Limb-wise widened a - b - c without normalization.
    BLST_FORCEINLINE void lufp_v2_sub2(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b, const UFpPack8x1& c)
    {
        for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
            r[limb] = VSUB(a[limb], b[limb]);
            r[limb] = VSUB(r[limb], c[limb]);
        }
    }

    // Forms a + b - p in widened radix-48 limbs before carry propagation.
    BLST_FORCEINLINE void v2_add_bias_neg_p(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b)
    {
        for (std::size_t i = 0; i < NWORDS; ++i) {
            r[i] = VADD(a[i], b[i]);
        }
        for (std::size_t i = 0; i < NWORDS; ++i) {
            r[NWORDS + i] = VSUB(VADD(a[NWORDS + i], b[NWORDS + i]), VSET1(static_cast<i64>(P48[i])));
        }
    }

    // Forms a - p in widened radix-48 limbs before carry propagation.
    BLST_FORCEINLINE void v2_bias_neg_p(UFpPack8x1& r, const UFpPack8x1& a)
    {
        for (std::size_t i = 0; i < NWORDS; ++i) {
            r[i] = a[i];
        }
        for (std::size_t i = 0; i < NWORDS; ++i) {
            r[NWORDS + i] = VSUB(a[NWORDS + i], VSET1(static_cast<i64>(P48[i])));
        }
    }

    // Returns the lane mask of widened values whose carried sign is negative.
    BLST_FORCEINLINE auto v2_neg_mask(const cc512i* limbs) -> v2_lane_mask_t;

    // Unsigned-style carry pass that masks each limb to the radix and forwards the carry.
    BLST_FORCEINLINE void v2_carry_range(cc512i* limbs, const std::size_t begin, const std::size_t end)
    {
        const cc512i bmask = VSET1(BMASK);
        for (std::size_t i = begin; i < end; ++i) {
            limbs[i + 1] = VADD(limbs[i + 1], VSRA(limbs[i], BRADIX));
            limbs[i] = VAND(limbs[i], bmask);
        }
    }

    // Signed carry pass over a fixed limb range in radix-48 representation.
    template <std::size_t Begin, std::size_t End>
    BLST_FORCEINLINE void v2_signed_carry_range(cc512i* limbs)
    {
        static_assert(Begin < (2 * NWORDS), "carry begin out of range");
        static_assert(End < (2 * NWORDS), "carry end out of range");
        static_assert(Begin <= End, "carry range must be non-decreasing");
        for (std::size_t i = Begin; i < End; ++i) {
            const cc512i carry = VSRA(limbs[i], BRADIX);
            limbs[i] = VSUB(limbs[i], VSHL(carry, BRADIX));
            limbs[i + 1] = VADD(limbs[i + 1], carry);
        }
    }

    // Carries a fixed limb range and returns the sign mask of the resulting widened values.
    template <std::size_t Begin, std::size_t End>
    BLST_FORCEINLINE auto v2_carry_range_with_sign(cc512i* limbs) -> v2_lane_mask_t
    {
        static_assert(Begin < (2 * NWORDS), "carry begin out of range");
        static_assert(End < (2 * NWORDS), "carry end out of range");
        static_assert(Begin <= End, "carry range must be non-decreasing");
        v2_carry_range(limbs, Begin, End);
        return v2_neg_mask(limbs);
    }

    // Legacy one-shot carry helper used by the older biased renorm experiments.
    BLST_FORCEINLINE auto v2_carry_from_biased(UFpPack8x1& r) -> cc512i
    {
        const cc512i bmask = VSET1(BMASK);

        r[1] = VADD(r[1], VSRA(r[0], BRADIX)); r[0] = VAND(r[0], bmask);
        r[2] = VADD(r[2], VSRA(r[1], BRADIX)); r[1] = VAND(r[1], bmask);
        r[3] = VADD(r[3], VSRA(r[2], BRADIX)); r[2] = VAND(r[2], bmask);
        r[4] = VADD(r[4], VSRA(r[3], BRADIX)); r[3] = VAND(r[3], bmask);
        r[5] = VADD(r[5], VSRA(r[4], BRADIX)); r[4] = VAND(r[4], bmask);
        r[6] = VADD(r[6], VSRA(r[5], BRADIX)); r[5] = VAND(r[5], bmask);
        r[7] = VADD(r[7], VSRA(r[6], BRADIX)); r[6] = VAND(r[6], bmask);
        r[8] = VADD(r[8], VSRA(r[7], BRADIX)); r[7] = VAND(r[7], bmask);
        cc512i t0 = VADD(r[9], VSRA(r[8], BRADIX));
        t0 = VADD(r[10], VSRA(t0, BRADIX));
        t0 = VADD(r[11], VSRA(t0, BRADIX));
        t0 = VADD(r[12], VSRA(t0, BRADIX));
        t0 = VADD(r[13], VSRA(t0, BRADIX));
        t0 = VADD(r[14], VSRA(t0, BRADIX));
        t0 = VADD(r[15], VSRA(t0, BRADIX));
        return VSRA(t0, 63);
    }

    // Adds one p on the negative lanes selected by smask and immediately re-carries the high half.
    BLST_FORCEINLINE void v2_repair_one_step_from_sign(UFpPack8x1& r, const cc512i smask)
    {
        const cc512i p0 = VSET1(P48[0]), p1 = VSET1(P48[1]), p2 = VSET1(P48[2]);
        const cc512i p3 = VSET1(P48[3]), p4 = VSET1(P48[4]), p5 = VSET1(P48[5]);
        const cc512i p6 = VSET1(P48[6]), p7 = VSET1(P48[7]);
        const cc512i bmask = VSET1(BMASK);

        r[8] = VADD(r[8], VAND(p0, smask)); r[9] = VADD(r[9], VAND(p1, smask));
        r[10] = VADD(r[10], VAND(p2, smask)); r[11] = VADD(r[11], VAND(p3, smask));
        r[12] = VADD(r[12], VAND(p4, smask)); r[13] = VADD(r[13], VAND(p5, smask));
        r[14] = VADD(r[14], VAND(p6, smask)); r[15] = VADD(r[15], VAND(p7, smask));

        r[9] = VADD(r[9], VSRA(r[8], BRADIX)); r[8] = VAND(r[8], bmask);
        r[10] = VADD(r[10], VSRA(r[9], BRADIX)); r[9] = VAND(r[9], bmask);
        r[11] = VADD(r[11], VSRA(r[10], BRADIX)); r[10] = VAND(r[10], bmask);
        r[12] = VADD(r[12], VSRA(r[11], BRADIX)); r[11] = VAND(r[11], bmask);
        r[13] = VADD(r[13], VSRA(r[12], BRADIX)); r[12] = VAND(r[12], bmask);
        r[14] = VADD(r[14], VSRA(r[13], BRADIX)); r[13] = VAND(r[13], bmask);
        r[15] = VADD(r[15], VSRA(r[14], BRADIX)); r[14] = VAND(r[14], bmask);
        r[15] = VAND(r[15], bmask);
    }

    // Generic sign-driven repair step: add Scale * p on the masked lanes and re-carry the high half.
    template <u64 Scale>
    BLST_FORCEINLINE auto v2_repair_sign_step(UFpPack8x1& r, const cc512i smask) -> cc512i
    {
        const auto& kP = v2_scaled_words<Scale>();
        const cc512i p0 = VSET1(kP[0]), p1 = VSET1(kP[1]), p2 = VSET1(kP[2]);
        const cc512i p3 = VSET1(kP[3]), p4 = VSET1(kP[4]), p5 = VSET1(kP[5]);
        const cc512i p6 = VSET1(kP[6]), p7 = VSET1(kP[7]);
        const cc512i bmask = VSET1(BMASK);

        r[8] = VADD(r[8], VAND(p0, smask)); r[9] = VADD(r[9], VAND(p1, smask));
        r[10] = VADD(r[10], VAND(p2, smask)); r[11] = VADD(r[11], VAND(p3, smask));
        r[12] = VADD(r[12], VAND(p4, smask)); r[13] = VADD(r[13], VAND(p5, smask));
        r[14] = VADD(r[14], VAND(p6, smask)); r[15] = VADD(r[15], VAND(p7, smask));

        r[9] = VADD(r[9], VSRA(r[8], BRADIX)); r[8] = VAND(r[8], bmask);
        r[10] = VADD(r[10], VSRA(r[9], BRADIX)); r[9] = VAND(r[9], bmask);
        r[11] = VADD(r[11], VSRA(r[10], BRADIX)); r[10] = VAND(r[10], bmask);
        r[12] = VADD(r[12], VSRA(r[11], BRADIX)); r[11] = VAND(r[11], bmask);
        r[13] = VADD(r[13], VSRA(r[12], BRADIX)); r[12] = VAND(r[12], bmask);
        r[14] = VADD(r[14], VSRA(r[13], BRADIX)); r[13] = VAND(r[13], bmask);
        r[15] = VADD(r[15], VSRA(r[14], BRADIX)); r[14] = VAND(r[14], bmask);
        r[15] = VAND(r[15], bmask);

        return VSRA(r[15], 63);
    }

    // Subtracts a fixed widened constant from the high half only.
    BLST_FORCEINLINE void v2_high_sub_const(cc512i* limbs, const std::array<cc512i, NWORDS>& words)
    {
        limbs[NWORDS + 0] = VSUB(limbs[NWORDS + 0], words[0]);
        limbs[NWORDS + 1] = VSUB(limbs[NWORDS + 1], words[1]);
        limbs[NWORDS + 2] = VSUB(limbs[NWORDS + 2], words[2]);
        limbs[NWORDS + 3] = VSUB(limbs[NWORDS + 3], words[3]);
        limbs[NWORDS + 4] = VSUB(limbs[NWORDS + 4], words[4]);
        limbs[NWORDS + 5] = VSUB(limbs[NWORDS + 5], words[5]);
        limbs[NWORDS + 6] = VSUB(limbs[NWORDS + 6], words[6]);
        limbs[NWORDS + 7] = VSUB(limbs[NWORDS + 7], words[7]);
    }

    // Adds a fixed widened constant to every lane of the high half.
    BLST_FORCEINLINE void v2_high_add_const(cc512i* limbs, const std::array<cc512i, NWORDS>& words)
    {
        limbs[NWORDS + 0] = VADD(limbs[NWORDS + 0], words[0]);
        limbs[NWORDS + 1] = VADD(limbs[NWORDS + 1], words[1]);
        limbs[NWORDS + 2] = VADD(limbs[NWORDS + 2], words[2]);
        limbs[NWORDS + 3] = VADD(limbs[NWORDS + 3], words[3]);
        limbs[NWORDS + 4] = VADD(limbs[NWORDS + 4], words[4]);
        limbs[NWORDS + 5] = VADD(limbs[NWORDS + 5], words[5]);
        limbs[NWORDS + 6] = VADD(limbs[NWORDS + 6], words[6]);
        limbs[NWORDS + 7] = VADD(limbs[NWORDS + 7], words[7]);
    }

    // Adds a fixed widened constant to the high half on the selected lanes.
    BLST_FORCEINLINE void v2_high_add_const_masked(
        cc512i* limbs, const std::array<cc512i, NWORDS>& words, const v2_lane_mask_t lane_mask)
    {
        limbs[NWORDS + 0] = VMADD(limbs[NWORDS + 0], lane_mask, limbs[NWORDS + 0], words[0]);
        limbs[NWORDS + 1] = VMADD(limbs[NWORDS + 1], lane_mask, limbs[NWORDS + 1], words[1]);
        limbs[NWORDS + 2] = VMADD(limbs[NWORDS + 2], lane_mask, limbs[NWORDS + 2], words[2]);
        limbs[NWORDS + 3] = VMADD(limbs[NWORDS + 3], lane_mask, limbs[NWORDS + 3], words[3]);
        limbs[NWORDS + 4] = VMADD(limbs[NWORDS + 4], lane_mask, limbs[NWORDS + 4], words[4]);
        limbs[NWORDS + 5] = VMADD(limbs[NWORDS + 5], lane_mask, limbs[NWORDS + 5], words[5]);
        limbs[NWORDS + 6] = VMADD(limbs[NWORDS + 6], lane_mask, limbs[NWORDS + 6], words[6]);
        limbs[NWORDS + 7] = VMADD(limbs[NWORDS + 7], lane_mask, limbs[NWORDS + 7], words[7]);
    }

    // Subtracts a fixed widened constant from the high half on the selected lanes.
    BLST_FORCEINLINE void v2_high_sub_const_masked(
        cc512i* limbs, const std::array<cc512i, NWORDS>& words, const v2_lane_mask_t lane_mask)
    {
        limbs[NWORDS + 0] = VMSUB(limbs[NWORDS + 0], lane_mask, limbs[NWORDS + 0], words[0]);
        limbs[NWORDS + 1] = VMSUB(limbs[NWORDS + 1], lane_mask, limbs[NWORDS + 1], words[1]);
        limbs[NWORDS + 2] = VMSUB(limbs[NWORDS + 2], lane_mask, limbs[NWORDS + 2], words[2]);
        limbs[NWORDS + 3] = VMSUB(limbs[NWORDS + 3], lane_mask, limbs[NWORDS + 3], words[3]);
        limbs[NWORDS + 4] = VMSUB(limbs[NWORDS + 4], lane_mask, limbs[NWORDS + 4], words[4]);
        limbs[NWORDS + 5] = VMSUB(limbs[NWORDS + 5], lane_mask, limbs[NWORDS + 5], words[5]);
        limbs[NWORDS + 6] = VMSUB(limbs[NWORDS + 6], lane_mask, limbs[NWORDS + 6], words[6]);
        limbs[NWORDS + 7] = VMSUB(limbs[NWORDS + 7], lane_mask, limbs[NWORDS + 7], words[7]);
    }

    // Lexicographically tests whether the carried high half is >= the provided constant.
    BLST_FORCEINLINE auto v2_high_ge_const_mask(
        const cc512i* limbs, const std::array<cc512i, NWORDS>& words) -> v2_lane_mask_t
    {
        v2_lane_mask_t eq_mask = static_cast<v2_lane_mask_t>(0xFF);
        v2_lane_mask_t gt_mask = static_cast<v2_lane_mask_t>(0);

        for (int i = static_cast<int>(NWORDS) - 1; i >= 0; --i) {
            const v2_lane_mask_t gt_i = v2_mask_gt_epi64(limbs[NWORDS + i], words[static_cast<std::size_t>(i)]);
            const v2_lane_mask_t eq_i = v2_mask_eq_epi64(limbs[NWORDS + i], words[static_cast<std::size_t>(i)]);
            gt_mask = static_cast<v2_lane_mask_t>(gt_mask | (eq_mask & gt_i));
            eq_mask = static_cast<v2_lane_mask_t>(eq_mask & eq_i);
        }

        return static_cast<v2_lane_mask_t>(gt_mask | eq_mask);
    }

    // u8 wrapper for the >= p test on the carried high half.
    BLST_FORCEINLINE auto v2_high_ge_p_mask_u8(const cc512i* limbs) -> u8
    {
        u8 eq_mask = 0xFF;
        u8 gt_mask = 0;

        for (int i = static_cast<int>(NWORDS) - 1; i >= 0; --i) {
            const cc512i mod = VSET1(static_cast<i64>(P48[static_cast<std::size_t>(i)]));
            const u8 gt_i = v2_mask_gt_epi64_u8(limbs[NWORDS + i], mod);
            const u8 eq_i = v2_mask_eq_epi64_u8(limbs[NWORDS + i], mod);
            gt_mask |= static_cast<u8>(eq_mask & gt_i);
            eq_mask &= eq_i;
        }

        return static_cast<u8>(gt_mask | eq_mask);
    }

    // u8 wrapper for the >= mod test on the carried high half.
    BLST_FORCEINLINE auto v2_high_ge_mod_mask_u8(const cc512i* limbs, const u64* mod) -> u8
    {
        u8 eq_mask = 0xFF;
        u8 gt_mask = 0;

        for (int i = static_cast<int>(NWORDS) - 1; i >= 0; --i) {
            const cc512i mod_i = VSET1(static_cast<i64>(mod[static_cast<std::size_t>(i)]));
            const u8 gt_i = v2_mask_gt_epi64_u8(limbs[NWORDS + i], mod_i);
            const u8 eq_i = v2_mask_eq_epi64_u8(limbs[NWORDS + i], mod_i);
            gt_mask |= static_cast<u8>(eq_mask & gt_i);
            eq_mask &= eq_i;
        }

        return static_cast<u8>(gt_mask | eq_mask);
    }

    // Adds a scalar radix-48 modulus on the selected lanes of the high half.
    BLST_FORCEINLINE void v2_high_add_mod_masked_u8(cc512i* limbs, const u8 lane_mask, const u64* mod)
    {
        for (std::size_t i = 0; i < NWORDS; ++i) {
            limbs[NWORDS + i] = VMADD(
                limbs[NWORDS + i], lane_mask, limbs[NWORDS + i], VSET1(static_cast<i64>(mod[i])));
        }
    }

    // Subtracts a scalar radix-48 modulus on the selected lanes of the high half.
    BLST_FORCEINLINE void v2_high_sub_mod_masked_u8(cc512i* limbs, const u8 lane_mask, const u64* mod)
    {
        for (std::size_t i = 0; i < NWORDS; ++i) {
            limbs[NWORDS + i] = VMSUB(
                limbs[NWORDS + i], lane_mask, limbs[NWORDS + i], VSET1(static_cast<i64>(mod[i])));
        }
    }

    // Reads the sign of the carried high half from the top limb.
    BLST_FORCEINLINE auto v2_neg_mask(const cc512i* limbs) -> v2_lane_mask_t
    {
        return v2_signbit_mask(limbs[(2 * NWORDS) - 1]);
    }

    // Exact widened equality, used only by debug validators.
    BLST_FORCEINLINE auto lufp_v2_equal(const UFpPack8x1& a, const UFpPack8x1& b) -> bool
    {
        for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
            if (v2_mask_eq_epi64_u8(a[limb], b[limb]) != 0xFFu) {
                return false;
            }
        }
        return true;
    }

    // Debug-only canonicalizer declarations used to state postconditions on fast renorm code.
    BLST_FORCEINLINE void renorm_any(UFpPack8x1& in_out);
    BLST_FORCEINLINE void renorm_any(UFp2Pack8x1x1& in_out);
    BLST_FORCEINLINE void renorm_any(UFp6Pack8x1& in_out);
    BLST_FORCEINLINE auto renorm_any_stable(const UFpPack8x1& got) -> bool;

    // Largest power of two <= value; used to choose the first strict-negative repair scale.
    BLST_FORCEINLINE constexpr auto v2_bit_floor_u32(const u32 value) -> u32
    {
        u32 floor = 0u;
        u32 x = value;
        while (x != 0u) {
            floor = (floor == 0u) ? 1u : (floor << 1u);
            x >>= 1u;
        }
        return floor;
    }

    // Fast sign-ladder primitive: add Scale * p on lanes that are already known negative.
    template <u64 Scale>
    BLST_FORCEINLINE auto v2_sign_step(cc512i* limbs, const v2_lane_mask_t neg_mask) -> v2_lane_mask_t
    {
        v2_high_add_const_masked(limbs, v2_scaled_vecs<Scale>(), neg_mask);
        return v2_carry_range_with_sign<NWORDS, (2 * NWORDS) - 1>(limbs);
    }

    // Runs a sequence of sign-driven repair steps, feeding the returned sign mask into the next step.
    template <u64... Scales>
    BLST_FORCEINLINE auto v2_sign_step_sequence(cc512i* limbs, v2_lane_mask_t sign_mask) -> v2_lane_mask_t
    {
        if constexpr (sizeof...(Scales) > 0) {
            ((sign_mask = v2_sign_step<Scales>(limbs, sign_mask)), ...);
        }
        return sign_mask;
    }

    // Subtracts an exact precomputed Scale * p from the high half.
    template <u32 Scale>
    BLST_FORCEINLINE void v2_high_sub_scaled(cc512i* limbs)
    {
        if constexpr (Scale == 0u) {
            return;
        }
        else {
            static_assert(Scale <= 16u, "scaled subtract requires direct precomputed p multiple");
            v2_high_sub_const(limbs, v2_scaled_vecs<Scale>());
        }
    }

    // Adds an exact precomputed Scale * p to every lane of the high half.
    template <u32 Scale>
    BLST_FORCEINLINE void v2_high_add_scaled(cc512i* limbs)
    {
        if constexpr (Scale == 0u) {
            return;
        }
        else {
            static_assert(Scale <= 16u, "scaled add requires direct precomputed p multiple");
            v2_high_add_const(limbs, v2_scaled_vecs<Scale>());
        }
    }

    // Adds an exact precomputed Scale * p to the selected lanes of the high half.
    template <u32 Scale>
    BLST_FORCEINLINE void v2_high_add_scaled_masked(cc512i* limbs, const v2_lane_mask_t lane_mask)
    {
        if (lane_mask == static_cast<v2_lane_mask_t>(0u)) {
            return;
        }
        if constexpr (Scale == 0u) {
            return;
        }
        else {
            static_assert(Scale <= 16u, "masked scaled add requires direct precomputed p multiple");
            v2_high_add_const_masked(limbs, v2_scaled_vecs<Scale>(), lane_mask);
        }
    }

    // Subtracts an exact precomputed Scale * p from the selected lanes of the high half.
    template <u32 Scale>
    BLST_FORCEINLINE void v2_high_sub_scaled_masked(cc512i* limbs, const v2_lane_mask_t lane_mask)
    {
        if (lane_mask == static_cast<v2_lane_mask_t>(0u)) {
            return;
        }
        if constexpr (Scale == 0u) {
            return;
        }
        else {
            static_assert(Scale <= 16u, "masked scaled subtract requires direct precomputed p multiple");
            v2_high_sub_const_masked(limbs, v2_scaled_vecs<Scale>(), lane_mask);
        }
    }

    // Reference strict-negative step: keep x + Scale * p only if the carried candidate is still negative.
    template <u32 Scale>
    BLST_FORCEINLINE void v2_keep_negative_step(UFpPack8x1& in_out)
    {
        v2_high_add_scaled<Scale>(in_out.data());
        const v2_lane_mask_t keep_mask = v2_carry_range_with_sign<NWORDS, (2 * NWORDS) - 1>(in_out.data());
        const v2_lane_mask_t revert_mask = static_cast<v2_lane_mask_t>(static_cast<u8>(0xFFu) ^ static_cast<u8>(keep_mask));
        if (revert_mask != static_cast<v2_lane_mask_t>(0u)) {
            v2_high_sub_scaled_masked<Scale>(in_out.data(), revert_mask);
            v2_signed_carry_range<NWORDS, (2 * NWORDS) - 1>(in_out.data());
        }
    }

    // Reference strict-negative ladder driven by descending powers of two.
    template <u32 Scale>
    BLST_FORCEINLINE void v2_keep_negative_ladder(UFpPack8x1& in_out)
    {
        static_assert(Scale != 0u, "repair ladder scale must be non-zero");
        v2_keep_negative_step<Scale>(in_out);
        if constexpr (Scale > 1u) {
            v2_keep_negative_ladder<Scale / 2u>(in_out);
        }
    }

    // Fused strict-negative step:
    //   x <- x + Scale * p - (!sign) * PrevScale * p
    // in one limb pass over the high half, then one carry/sign pass.
    // PrevScale == 0 seeds the ladder without any conditional subtract.
    template <u32 PrevScale, u32 Scale>
    BLST_FORCEINLINE auto v2_keep_negative_fuse_step(
        cc512i* limbs, const v2_lane_mask_t sign) -> v2_lane_mask_t
    {
        const auto& scale_words = v2_scaled_vecs<Scale>();
        if constexpr (PrevScale == 0u) {
            for (std::size_t i = 0; i < NWORDS; ++i) {
                limbs[NWORDS + i] = VADD(limbs[NWORDS + i], scale_words[i]);
            }
        }
        else {
            const auto& prev_words = v2_scaled_vecs<PrevScale>();
            const v2_lane_mask_t nonsign =
                static_cast<v2_lane_mask_t>(static_cast<u8>(0xFFu) ^ static_cast<u8>(sign));
            for (std::size_t i = 0; i < NWORDS; ++i) {
                const cc512i candidate = VADD(limbs[NWORDS + i], scale_words[i]);
                limbs[NWORDS + i] = VMSUB(candidate, nonsign, candidate, prev_words[i]);
            }
        }
        return v2_carry_range_with_sign<NWORDS, (2 * NWORDS) - 1>(limbs);
    }

    // Recursive fused ladder. The base case <1,1> is the folded final +p canonicalization step.
    template <u32 PrevScale, u32 Scale>
    BLST_FORCEINLINE auto v2_keep_negative_fuse_ladder(
        cc512i* limbs, const v2_lane_mask_t sign) -> v2_lane_mask_t
    {
        const v2_lane_mask_t next_sign = v2_keep_negative_fuse_step<PrevScale, Scale>(limbs, sign);
        if constexpr (Scale > 1u) {
            return v2_keep_negative_fuse_ladder<Scale, Scale / 2u>(limbs, next_sign);
        }
        else {
            return v2_keep_negative_fuse_step<1u, 1u>(limbs, next_sign);
        }
    }

    // Proven strict-negative renorm:
    //   bias by B * p, repair with keep-negative steps, then add the final +p explicitly.
    template <u32 A, u32 B>
    BLST_FORCEINLINE void lufp_v2_renorm_negative_fold(UFpPack8x1& in_out)
    {
        static_assert((A + B) > 0u, "repair source class must have non-zero width");

        constexpr u32 kBias = B;
        constexpr u32 kFirstScale = v2_bit_floor_u32(A + B);

        // Bias into a strict-negative class, then keep only still-negative repaired candidates.
        v2_high_sub_scaled<kBias>(in_out.data());
        v2_signed_carry_range<0, (2 * NWORDS) - 1>(in_out.data());
        if constexpr (kFirstScale != 0u) {
            v2_keep_negative_ladder<kFirstScale>(in_out);
        }

        // Unbias into the final non-negative canonical class.
        v2_high_add_scaled<1u>(in_out.data());
        v2_signed_carry_range<NWORDS, (2 * NWORDS) - 1>(in_out.data());
        BLST_ASSERT(renorm_any_stable(in_out));
    }

    // Fused strict-negative renorm for a proved one-sided source class x in [-A*p, B*p).
    //
    // Reference algorithm:
    //   1. bias into strict negative: x <- x - B*p, so x in [-(A+B)*p, 0)
    //   2. walk descending repair scales a_i = bit_floor(A+B), ..., 2, 1
    //      keeping only candidates that remain negative:
    //          x' = x + a_i*p
    //          if x' < 0: x <- x'
    //   3. fold the final +p canonicalization into the terminal <1,1> step
    //
    // The fused machine implements the same state transition, but carries forward the
    // previous keep/reject decision as a sign bit and updates the next candidate by:
    //     x <- x + a_i*p - (!sign) * a_{i-1}*p
    // so a rejected prior candidate is undone implicitly instead of by restoring an
    // older state.
    //
    // The function now has two entry modes:
    //   - seeded fast path, used when B <= first_scale:
    //       instead of materializing x - B*p and then applying the first step +first_scale*p,
    //       start directly from their sum:
    //           (x - B*p) + first_scale*p = x + (first_scale - B)*p
    //       After the initial carry/sign pass, the state is exactly the post-first-step
    //       candidate of the reference algorithm, so the remaining fused ladder can start
    //       at PrevScale = first_scale.
    //   - generic fallback:
    //       perform the explicit bias x <- x - B*p first, then run the original fused ladder
    //       seeded from PrevScale = 0.
    //
    // This is still the same strict-negative proof contract as
    // lufp_v2_renorm_negative_fold(...); only the way the first bias/step pair is
    // materialized has changed.
    template <u32 A, u32 B>
    BLST_FORCEINLINE void lufp_v2_renorm_negative_fold_fused(UFpPack8x1& in_out)
    {
        static_assert((A + B) > 0u, "repair source class must have non-zero width");

        constexpr u32 kBias = B;
        constexpr u32 kFirstScale = v2_bit_floor_u32(A + B);

        if constexpr (kFirstScale != 0u) {
            if constexpr (kBias <= kFirstScale) {
                // Start directly at the post-first-step candidate:
                // (x - B*p) + kFirstScale*p == x + (kFirstScale - B)*p.
                constexpr u32 kSeedScale = kFirstScale - kBias;
                v2_high_add_scaled<kSeedScale>(in_out.data());
                const v2_lane_mask_t sign_mask =
                    v2_carry_range_with_sign<0, (2 * NWORDS) - 1>(in_out.data());
                if constexpr (kFirstScale > 1u) {
                    static_cast<void>(
                        v2_keep_negative_fuse_ladder<kFirstScale, kFirstScale / 2u>(in_out.data(), sign_mask));
                }
                else {
                    static_cast<void>(v2_keep_negative_fuse_step<1u, 1u>(in_out.data(), sign_mask));
                }
            }
            else {
                v2_high_sub_scaled<kBias>(in_out.data());
                v2_signed_carry_range<0, (2 * NWORDS) - 1>(in_out.data());
                static_cast<void>(
                    v2_keep_negative_fuse_ladder<0u, kFirstScale>(in_out.data(), static_cast<v2_lane_mask_t>(0u)));
            }
        }
        else {
            v2_high_sub_scaled<kBias>(in_out.data());
            v2_signed_carry_range<0, (2 * NWORDS) - 1>(in_out.data());
        }

        BLST_ASSERT(renorm_any_stable(in_out));
    }

    // Very slow canonicalizer used only as a debug oracle for widened half values.
    BLST_FORCEINLINE void renorm_any(UFpPack8x1& in_out)
    {
        const cc512i bmask = VSET1(BMASK);

        v2_signed_carry_range<0, (2 * NWORDS) - 1>(in_out.data());

        for (std::size_t iter = 0; iter < 256; ++iter) {
            const u8 neg_mask = static_cast<u8>(v2_signbit_mask(in_out[(2 * NWORDS) - 1]));
            if (neg_mask != 0u) {
                v2_high_add_mod_masked_u8(in_out.data(), neg_mask, kV2P48X4.data());
                v2_signed_carry_range<NWORDS, (2 * NWORDS) - 1>(in_out.data());
                continue;
            }

            const u8 ge_4p_mask = v2_high_ge_mod_mask_u8(in_out.data(), kV2P48X4.data());
            if (ge_4p_mask != 0u) {
                v2_high_sub_mod_masked_u8(in_out.data(), ge_4p_mask, kV2P48X4.data());
                v2_signed_carry_range<NWORDS, (2 * NWORDS) - 1>(in_out.data());
                continue;
            }

            const u8 ge_2p_mask = v2_high_ge_mod_mask_u8(in_out.data(), kV2P48X2.data());
            if (ge_2p_mask != 0u) {
                v2_high_sub_mod_masked_u8(in_out.data(), ge_2p_mask, kV2P48X2.data());
                v2_signed_carry_range<NWORDS, (2 * NWORDS) - 1>(in_out.data());
                continue;
            }

            const u8 ge_p_mask = v2_high_ge_p_mask_u8(in_out.data());
            if (ge_p_mask != 0u) {
                v2_high_sub_mod_masked_u8(in_out.data(), ge_p_mask, kV2P48X1.data());
                v2_signed_carry_range<NWORDS, (2 * NWORDS) - 1>(in_out.data());
                continue;
            }

            for (std::size_t limb = 0; limb < (2 * NWORDS); ++limb) {
                in_out[limb] = VAND(in_out[limb], bmask);
            }
            return;
        }

        std::fprintf(stderr, "renorm_any failed to converge\n");
        std::abort();
    }

    // Debug oracle lifted to Fp2.
    BLST_FORCEINLINE void renorm_any(UFp2Pack8x1x1& in_out)
    {
        renorm_any(in_out[0]);
        renorm_any(in_out[1]);
    }

    // Debug oracle lifted to Fp6.
    BLST_FORCEINLINE void renorm_any(UFp6Pack8x1& in_out)
    {
        renorm_any(in_out.c[0]);
        renorm_any(in_out.c[1]);
        renorm_any(in_out.c[2]);
    }

    // Checks whether a fast renorm result is already a fixed point of renorm_any.
    BLST_FORCEINLINE auto renorm_any_stable(const UFpPack8x1& got) -> bool
    {
        UFpPack8x1 canonical = got;
        renorm_any(canonical);
        return lufp_v2_equal(got, canonical);
    }

    // Single-step exact-add repair used by some historical relaxed renorm experiments.
    BLST_FORCEINLINE void lufp_v2_renorm_exact_add(UFpPack8x1& in_out)
    {
        const cc512i smask = v2_carry_from_biased(in_out);
        v2_repair_one_step_from_sign(in_out, smask);
    }

    // Historical common FP6 half renorm: carry, fix negatives, then trim back from [0, 4p) to [0, p).
    BLST_FORCEINLINE void lufp_v2_renorm_fp6_mul_common(UFpPack8x1& in_out)
    {
        const cc512i bmask = VSET1(BMASK);

        v2_carry_range(in_out.data(), 0, (2 * NWORDS) - 1);

        const v2_lane_mask_t neg_mask = v2_mask_gt_epi64(VZERO(), in_out[(2 * NWORDS) - 1]);
        v2_high_add_const_masked(in_out.data(), kV2P48X4Vec, neg_mask);
        v2_carry_range(in_out.data(), NWORDS, (2 * NWORDS) - 1);

        const v2_lane_mask_t neg_fix_mask = v2_mask_gt_epi64(VZERO(), in_out[(2 * NWORDS) - 1]);
        v2_high_add_const_masked(in_out.data(), kV2P48X1Vec, neg_fix_mask);
        v2_carry_range(in_out.data(), NWORDS, (2 * NWORDS) - 1);

        const v2_lane_mask_t ge_4p_mask = v2_high_ge_const_mask(in_out.data(), kV2P48X4Vec);
        v2_high_sub_const_masked(in_out.data(), kV2P48X4Vec, ge_4p_mask);
        v2_carry_range(in_out.data(), NWORDS, (2 * NWORDS) - 1);

        const v2_lane_mask_t ge_2p_mask = v2_high_ge_const_mask(in_out.data(), kV2P48X2Vec);
        v2_high_sub_const_masked(in_out.data(), kV2P48X2Vec, ge_2p_mask);
        v2_carry_range(in_out.data(), NWORDS, (2 * NWORDS) - 1);

        const v2_lane_mask_t ge_p_mask = v2_high_ge_const_mask(in_out.data(), kV2P48X1Vec);
        v2_high_sub_const_masked(in_out.data(), kV2P48X1Vec, ge_p_mask);
        v2_carry_range(in_out.data(), NWORDS, (2 * NWORDS) - 1);

        for (std::size_t limb = 0; limb < (2 * NWORDS); ++limb) {
            in_out[limb] = VAND(in_out[limb], bmask);
        }
    }

    //struct UFp2Pack8x1x1
    //{
    //    UFpPack8x1 c[2];
    //
    //    BLST_FORCEINLINE UFpPack8x1& operator[](std::size_t i) { return c[i]; }
    //    BLST_FORCEINLINE const UFpPack8x1& operator[](std::size_t i) const { return c[i]; }
    //};
    //
    //static_assert(sizeof(UFp2Pack8x1x1) == sizeof(UFp2Pack8x1x1));
    //static_assert(alignof(UFp2Pack8x1x1) == alignof(UFp2Pack8x1x1));

    //BLST_FORCEINLINE auto lufp_v2_alias(UFp2Pack8x1x1& src) -> UFp2Pack8x1x1&
    //{
    //    return reinterpret_cast<UFp2Pack8x1x1&>(src);
    //}
    //
    //BLST_FORCEINLINE auto lufp_v2_alias(const UFp2Pack8x1x1& src) -> const UFp2Pack8x1x1&
    //{
    //    return reinterpret_cast<const UFp2Pack8x1x1&>(src);
    //}

    // Widened Fp2 add done component-wise.
    BLST_FORCEINLINE void lufp2_v2_add(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
    {
        lufp_v2_add(r[0], a[0], b[0]);
        lufp_v2_add(r[1], a[1], b[1]);
    }

    // Widened Fp2 subtract done component-wise.
    BLST_FORCEINLINE void lufp2_v2_sub(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b)
    {
        lufp_v2_sub(r[0], a[0], b[0]);
        lufp_v2_sub(r[1], a[1], b[1]);
    }

    // Widened Fp2 a - b - c done component-wise.
    BLST_FORCEINLINE void lufp2_v2_sub2(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b, const UFp2Pack8x1x1& c)
    {
        lufp_v2_sub2(r[0], a[0], b[0], c[0]);
        lufp_v2_sub2(r[1], a[1], b[1], c[1]);
    }

    // Widened Fp2 multiply by (u + 1), used in the tower formulas.
    BLST_FORCEINLINE void lufp2_v2_mul_by_u_plus_1(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a)
    {
        lufp_v2_sub(r[0], a[0], a[1]);
        lufp_v2_add(r[1], a[0], a[1]);
    }

    // Widened Fp6 add done component-wise with no repair.
    BLST_FORCEINLINE void lufp6_v2_add(UFp6Pack8x1& r, const UFp6Pack8x1& a, const UFp6Pack8x1& b)
    {
        lufp2_v2_add(r.c[0], a.c[0], b.c[0]);
        lufp2_v2_add(r.c[1], a.c[1], b.c[1]);
        lufp2_v2_add(r.c[2], a.c[2], b.c[2]);
    }

    // Widened Fp6 subtract done component-wise with no repair.
    BLST_FORCEINLINE void lufp6_v2_sub(UFp6Pack8x1& r, const UFp6Pack8x1& a, const UFp6Pack8x1& b)
    {
        lufp2_v2_sub(r.c[0], a.c[0], b.c[0]);
        lufp2_v2_sub(r.c[1], a.c[1], b.c[1]);
        lufp2_v2_sub(r.c[2], a.c[2], b.c[2]);
    }

    // Widened Fp6 subtract of two inputs done component-wise with no repair.
    BLST_FORCEINLINE void lufp6_v2_sub2(UFp6Pack8x1& r, const UFp6Pack8x1& a, const UFp6Pack8x1& b, const UFp6Pack8x1& c)
    {
        lufp2_v2_sub2(r.c[0], a.c[0], b.c[0], c.c[0]);
        lufp2_v2_sub2(r.c[1], a.c[1], b.c[1], c.c[1]);
        lufp2_v2_sub2(r.c[2], a.c[2], b.c[2], c.c[2]);
    }

    // Widened Fp6 multiply by the tower nonresidue with no repair.
    BLST_FORCEINLINE void lufp6_v2_mul_by_nonresidue(UFp6Pack8x1& r, const UFp6Pack8x1& a)
    {
        lufp2_v2_mul_by_u_plus_1(r.c[0], a.c[2]);
        r.c[1] = a.c[0];
        r.c[2] = a.c[1];
    }

    // Widened Fp6 r = a + v*b with no repair.
    BLST_FORCEINLINE void lufp6_v2_add_mul_by_nonresidue(UFp6Pack8x1& r, const UFp6Pack8x1& a, const UFp6Pack8x1& b)
    {
        lufp2_v2_mul_by_u_plus_1(r.c[0], b.c[2]);
        lufp2_v2_add(r.c[0], a.c[0], r.c[0]);
        lufp2_v2_add(r.c[1], a.c[1], b.c[0]);
        lufp2_v2_add(r.c[2], a.c[2], b.c[1]);
    }

    // Carries only the high half and returns the resulting sign mask.
    BLST_FORCEINLINE auto v2_carry_high_half_with_sign(cc512i* limbs) -> v2_lane_mask_t
    {
        const cc512i bmask = VSET1(BMASK);
        for (std::size_t i = 0; i < NWORDS - 1; ++i) {
            limbs[i + 1] = VADD(limbs[i + 1], VSRA(limbs[i], BRADIX));
            limbs[i] = VAND(limbs[i], bmask);
        }
        limbs[NWORDS - 1] = VAND(limbs[NWORDS - 1], bmask);
        return v2_signbit_mask(limbs[NWORDS - 1]);
    }

    //// Historical fast FP6 half renorm used by the older common-path experiments.
    //BLST_FORCEINLINE void lufp_v2_fp6_mul_common_renorm_half_fast(UFpPack8x1& r)
    //{
    //    v2_high_sub_const(r.data(), kV2Bias3PPlus3Vec);
    //    v2_lane_mask_t sign_mask = v2_carry_range_with_sign<0, (2 * NWORDS) - 1>(r.data());
    //    sign_mask = v2_sign_step<4>(r.data(), sign_mask);
    //    sign_mask = v2_sign_step<2>(r.data(), sign_mask);
    //    sign_mask = v2_sign_step<1>(r.data(), sign_mask);
    //    v2_sign_step<1>(r.data(), sign_mask);
    //    BLST_ASSERT(renorm_any_stable(r));
    //}

    //// Historical step-candidate FP6 half renorm kept for debugging alternative carry/sign orders.
    //BLST_FORCEINLINE void lufp_v2_fp6_mul_common_renorm_half_step_candidate(UFpPack8x1& r)
    //{
    //    v2_high_sub_const(r.data(), kV2Bias3PPlus3Vec);
    //    v2_signed_carry_range<0, (2 * NWORDS) - 1>(r.data());
    //    v2_lane_mask_t sign_mask = v2_signbit_mask(r[(2 * NWORDS) - 1]);
    //    sign_mask = v2_sign_step<4>(r.data(), sign_mask);
    //    sign_mask = v2_sign_step<2>(r.data(), sign_mask);
    //    sign_mask = v2_sign_step<1>(r.data(), sign_mask);
    //    v2_sign_step<1>(r.data(), sign_mask);
    //    BLST_ASSERT(renorm_any_stable(r));
    //}

    //// Conservative FP6 half renorm used as a correctness oracle for the common path.
    //BLST_FORCEINLINE void lufp_v2_fp6_mul_common_renorm_half_safe(UFpPack8x1& r)
    //{
    //    const cc512i bmask = VSET1(BMASK);

    //    v2_signed_carry_range<0, (2 * NWORDS) - 1>(r.data());

    //    const u8 neg_mask = static_cast<u8>(v2_signbit_mask(r[(2 * NWORDS) - 1]));
    //    v2_high_add_mod_masked_u8(r.data(), neg_mask, kV2P48X4.data());
    //    v2_signed_carry_range<NWORDS, (2 * NWORDS) - 1>(r.data());

    //    const u8 neg_fix_mask = static_cast<u8>(v2_signbit_mask(r[(2 * NWORDS) - 1]));
    //    v2_high_add_mod_masked_u8(r.data(), neg_fix_mask, kV2P48X1.data());
    //    v2_signed_carry_range<NWORDS, (2 * NWORDS) - 1>(r.data());

    //    const u8 ge_4p_mask = v2_high_ge_mod_mask_u8(r.data(), kV2P48X4.data());
    //    v2_high_sub_mod_masked_u8(r.data(), ge_4p_mask, kV2P48X4.data());
    //    v2_signed_carry_range<NWORDS, (2 * NWORDS) - 1>(r.data());

    //    const u8 ge_2p_mask = v2_high_ge_mod_mask_u8(r.data(), kV2P48X2.data());
    //    v2_high_sub_mod_masked_u8(r.data(), ge_2p_mask, kV2P48X2.data());
    //    v2_signed_carry_range<NWORDS, (2 * NWORDS) - 1>(r.data());

    //    const u8 ge_p_mask = v2_high_ge_p_mask_u8(r.data());
    //    v2_high_sub_mod_masked_u8(r.data(), ge_p_mask, kV2P48X1.data());
    //    v2_signed_carry_range<NWORDS, (2 * NWORDS) - 1>(r.data());

    //    for (std::size_t limb = 0; limb < (2 * NWORDS); ++limb) {
    //        r[limb] = VAND(r[limb], bmask);
    //    }
    //}

    //// Historical c0 real-half combine + fast renorm helper for the common-path experiments.
    //BLST_FORCEINLINE void lufp_v2_c0_lo_renorm(
    //    UFpPack8x1& r,
    //    const UFpPack8x1& a0,
    //    const UFpPack8x1& b0,
    //    const UFpPack8x1& c0,
    //    const UFpPack8x1& a1,
    //    const UFpPack8x1& b1,
    //    const UFpPack8x1& c1,
    //    const UFpPack8x1& d0)
    //{
    //    for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
    //        const cc512i x0 = VSUB(VSUB(a0[limb], b0[limb]), c0[limb]);
    //        const cc512i x1 = VSUB(VSUB(a1[limb], b1[limb]), c1[limb]);
    //        r[limb] = VADD(VSUB(x0, x1), d0[limb]);
    //    }
    //    lufp_v2_fp6_mul_common_renorm_half_fast(r);
    //}

    //// Historical c0 imaginary-half combine + fast renorm helper for the common-path experiments.
    //BLST_FORCEINLINE void lufp_v2_c0_hi_renorm(
    //    UFpPack8x1& r,
    //    const UFpPack8x1& a0,
    //    const UFpPack8x1& b0,
    //    const UFpPack8x1& c0,
    //    const UFpPack8x1& a1,
    //    const UFpPack8x1& b1,
    //    const UFpPack8x1& c1,
    //    const UFpPack8x1& d1)
    //{
    //    for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
    //        const cc512i x0 = VSUB(VSUB(a0[limb], b0[limb]), c0[limb]);
    //        const cc512i x1 = VSUB(VSUB(a1[limb], b1[limb]), c1[limb]);
    //        r[limb] = VADD(VADD(x0, x1), d1[limb]);
    //    }
    //    lufp_v2_fp6_mul_common_renorm_half_fast(r);
    //}

    //// Historical c1 real-half combine + fast renorm helper for the common-path experiments.
    //BLST_FORCEINLINE void lufp_v2_c1_lo_renorm(
    //    UFpPack8x1& r,
    //    const UFpPack8x1& a,
    //    const UFpPack8x1& b,
    //    const UFpPack8x1& c,
    //    const UFpPack8x1& d0,
    //    const UFpPack8x1& d1)
    //{
    //    for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
    //        r[limb] = VADD(VSUB(VSUB(a[limb], b[limb]), c[limb]), VSUB(d0[limb], d1[limb]));
    //    }
    //    lufp_v2_fp6_mul_common_renorm_half_fast(r);
    //}

    //// Historical c1 imaginary-half combine + fast renorm helper for the common-path experiments.
    //BLST_FORCEINLINE void lufp_v2_c1_hi_renorm(
    //    UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b, const UFpPack8x1& c, const UFpPack8x1& d0, const UFpPack8x1& d1)
    //{
    //    for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
    //        r[limb] = VADD(VSUB(VSUB(a[limb], b[limb]), c[limb]), VADD(d0[limb], d1[limb]));
    //    }
    //    lufp_v2_fp6_mul_common_renorm_half_fast(r);
    //}

    //// Historical fused subtract/add/fast-renorm helper for the common-path experiments.
    //BLST_FORCEINLINE void lufp_v2_sub2_add_renorm(
    //    UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b, const UFpPack8x1& c, const UFpPack8x1& d)
    //{
    //    for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
    //        r[limb] = VADD(VSUB(VSUB(a[limb], b[limb]), c[limb]), d[limb]);
    //    }
    //    lufp_v2_fp6_mul_common_renorm_half_fast(r);
    //}

    //// Applies the historical fast common FP6 renorm to both Fp halves.
    //BLST_FORCEINLINE void lufp2_v2_renorm_fp6_mul_common_fast(UFp2Pack8x1x1& in_out)
    //{
    //    lufp_v2_fp6_mul_common_renorm_half_fast(in_out[0]);
    //    lufp_v2_fp6_mul_common_renorm_half_fast(in_out[1]);
    //}

    //// Applies the historical final-step candidate common FP6 renorm to both Fp halves.
    //BLST_FORCEINLINE void lufp2_v2_renorm_fp6_mul_common_final_step(UFp2Pack8x1x1& in_out)
    //{
    //    lufp_v2_fp6_mul_common_renorm_half_step_candidate(in_out[0]);
    //    lufp_v2_fp6_mul_common_renorm_half_step_candidate(in_out[1]);
    //}

    //// Applies the historical step-candidate common FP6 renorm to both Fp halves.
    //BLST_FORCEINLINE void lufp2_v2_renorm_fp6_mul_common_step_candidate(UFp2Pack8x1x1& in_out)
    //{
    //    lufp_v2_fp6_mul_common_renorm_half_step_candidate(in_out[0]);
    //    lufp_v2_fp6_mul_common_renorm_half_step_candidate(in_out[1]);
    //}

    //// Applies the conservative common FP6 renorm to both Fp halves.
    //BLST_FORCEINLINE void lufp2_v2_renorm_fp6_mul_common_safe(UFp2Pack8x1x1& in_out)
    //{
    //    lufp_v2_fp6_mul_common_renorm_half_safe(in_out[0]);
    //    lufp_v2_fp6_mul_common_renorm_half_safe(in_out[1]);
    //}

    //// Applies the exact-add repair to both Fp halves.
    //BLST_FORCEINLINE void lufp2_v2_renorm_exact_add(UFp2Pack8x1x1& in_out)
    //{
    //    lufp_v2_renorm_exact_add(in_out[0]);
    //    lufp_v2_renorm_exact_add(in_out[1]);
    //}


    template <bool DoRenorm = true>
    inline void mul_fp6_pack_lufpv2_wide(UFp6Pack8x1& r, const Fp6Pack8x1& a, const Fp6Pack8x1& b);
    inline void mul_fp12_pack_lufpv2_wide(UFp12Pack8x1& r, const Fp12Pack8x1& a, const Fp12Pack8x1& b);

    inline void mul_fp12_pack_lufpv2_correct(UFp12Pack8x1& r, const Fp12Pack8x1& a, const Fp12Pack8x1& b)
    {
        mul_fp12_pack_lufpv2_wide(r, a, b);
    }
    // Baseline FP6 LUFp-v2 multiply: raw combine followed by the fused strict-negative renorm.
    template <bool DoRenorm>
    inline void mul_fp6_pack_lufpv2_wide(UFp6Pack8x1& r, const Fp6Pack8x1& a, const Fp6Pack8x1& b)
    {
        Fp2Pack8x1x1 aa, bb;
        UFp2Pack8x1x1 t0, t1, t2;
        UFp2Pack8x1x1 w;

        UFp2Pack8x1x1::mul(t0, a.c[0], b.c[0]);
        UFp2Pack8x1x1::mul(t1, a.c[1], b.c[1]);
        UFp2Pack8x1x1::mul(t2, a.c[2], b.c[2]);

        Fp2Pack8x1x1::add(aa, a.c[1], a.c[2]);
        Fp2Pack8x1x1::add(bb, b.c[1], b.c[2]);
        UFp2Pack8x1x1::mul(w, aa, bb);
        for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
            const cc512i x00 = VSUB(VSUB(w[0][limb], t1[0][limb]), t2[0][limb]);
            const cc512i x01 = VSUB(VSUB(w[1][limb], t1[1][limb]), t2[1][limb]);
            r.c[0][0][limb] = VADD(VSUB(x00, x01), t0[0][limb]);

            const cc512i x10 = VSUB(VSUB(w[0][limb], t1[0][limb]), t2[0][limb]);
            const cc512i x11 = VSUB(VSUB(w[1][limb], t1[1][limb]), t2[1][limb]);
            r.c[0][1][limb] = VADD(VADD(x10, x11), t0[1][limb]);
        }
        if constexpr (DoRenorm) {
            lufp_v2_renorm_negative_fold_fused<3u, 4u>(r.c[0][0]);
            lufp_v2_renorm_negative_fold_fused<4u, 3u>(r.c[0][1]);
        }

        Fp2Pack8x1x1::add(aa, a.c[0], a.c[1]);
        Fp2Pack8x1x1::add(bb, b.c[0], b.c[1]);
        UFp2Pack8x1x1::mul(w, aa, bb);
        for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
            r.c[1][0][limb] = VADD(
                VSUB(VSUB(w[0][limb], t0[0][limb]), t1[0][limb]),
                VSUB(t2[0][limb], t2[1][limb]));
            r.c[1][1][limb] = VADD(
                VSUB(VSUB(w[1][limb], t0[1][limb]), t1[1][limb]),
                VADD(t2[0][limb], t2[1][limb]));
        }
        if constexpr (DoRenorm) {
            lufp_v2_renorm_negative_fold_fused<3u, 2u>(r.c[1][0]);
            lufp_v2_renorm_negative_fold_fused<2u, 3u>(r.c[1][1]);
        }

        Fp2Pack8x1x1::add(aa, a.c[0], a.c[2]);
        Fp2Pack8x1x1::add(bb, b.c[0], b.c[2]);
        UFp2Pack8x1x1::mul(w, aa, bb);
        for (std::size_t limb = 0; limb < 2 * NWORDS; ++limb) {
            r.c[2][0][limb] = VADD(VSUB(VSUB(w[0][limb], t0[0][limb]), t2[0][limb]), t1[0][limb]);
            r.c[2][1][limb] = VADD(VSUB(VSUB(w[1][limb], t0[1][limb]), t2[1][limb]), t1[1][limb]);
        }
        if constexpr (DoRenorm) {
            lufp_v2_renorm_negative_fold_fused<2u, 2u>(r.c[2][0]);
            lufp_v2_renorm_negative_fold_fused<2u, 2u>(r.c[2][1]);
        }
    }

    // True always-wide Fp12 candidate: no carries or renorms after the inner Fp6 multiplies.
    inline void mul_fp12_pack_lufpv2_wide(UFp12Pack8x1& r, const Fp12Pack8x1& a, const Fp12Pack8x1& b)
    {
        UFp6Pack8x1 t0, t1;

        mul_fp6_pack_lufpv2_wide<false>(t0, a.c[0], b.c[0]);
        mul_fp6_pack_lufpv2_wide<false>(t1, a.c[1], b.c[1]);

        {
            Fp6Pack8x1 t3, t4;
            Fp6Pack8x1::add(t3, a.c[0], a.c[1]);
            Fp6Pack8x1::add(t4, b.c[0], b.c[1]);
            mul_fp6_pack_lufpv2_wide<false>(r.c[1], t3, t4);
        }

        lufp6_v2_sub2(r.c[1], r.c[1], t0, t1);
        lufp_v2_renorm_negative_fold_fused<11u, 10u>(r.c[1].c[0][0]);
        lufp_v2_renorm_negative_fold_fused<10u, 11u>(r.c[1].c[0][1]);
        lufp_v2_renorm_negative_fold_fused<7u, 8u>(r.c[1].c[1][0]);
        lufp_v2_renorm_negative_fold_fused<8u, 7u>(r.c[1].c[1][1]);
        lufp_v2_renorm_negative_fold_fused<6u, 6u>(r.c[1].c[2][0]);
        lufp_v2_renorm_negative_fold_fused<6u, 6u>(r.c[1].c[2][1]);

        lufp6_v2_add_mul_by_nonresidue(r.c[0], t0, t1);

        // Final Fp12 per-half source classes from the worksheet.
        lufp_v2_renorm_negative_fold_fused<7u, 8u>(r.c[0].c[0][0]);
        lufp_v2_renorm_negative_fold_fused<8u, 7u>(r.c[0].c[0][1]);
        lufp_v2_renorm_negative_fold_fused<6u, 6u>(r.c[0].c[1][0]);
        lufp_v2_renorm_negative_fold_fused<6u, 6u>(r.c[0].c[1][1]);
        lufp_v2_renorm_negative_fold_fused<5u, 4u>(r.c[0].c[2][0]);
        lufp_v2_renorm_negative_fold_fused<4u, 5u>(r.c[0].c[2][1]);
    }

    inline void mul_fp12_pack_lufpv2_fp6_safe_wide(UFp12Pack8x1& r, const Fp12Pack8x1& a, const Fp12Pack8x1& b)
    {
        UFp6Pack8x1 t0, t1, t2;
        Fp6Pack8x1 t3, t4;

        mul_fp6_pack_lufpv2_wide<true>(t0, a.c[0], b.c[0]);
        mul_fp6_pack_lufpv2_wide<true>(t1, a.c[1], b.c[1]);

        Fp6Pack8x1::add(t3, a.c[0], a.c[1]);
        Fp6Pack8x1::add(t4, b.c[0], b.c[1]);
        mul_fp6_pack_lufpv2_wide<true>(r.c[1], t3, t4);
        UFp6Pack8x1::sub(r.c[1], r.c[1], t0);
        UFp6Pack8x1::sub(r.c[1], r.c[1], t1);

        UFp6Pack8x1::mul_by_nonresidue(t2, t1);
        UFp6Pack8x1::add(r.c[0], t0, t2);
    }

}
