#!/usr/bin/env python3
# Copyright 2026 Category Labs, Inc.
# SPDX-License-Identifier: Apache-2.0

"""Fail closed unless production AVX dispatch selects the proved MatrixNoK owners."""

from __future__ import annotations

from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def compact(path: Path) -> str:
    require(path.is_file(), f"missing dispatch owner: {path}")
    return " ".join(path.read_text(encoding="utf-8").split())


def require_fragments(path: Path, fragments: tuple[str, ...]) -> None:
    source = compact(path)
    for fragment in fragments:
        require(fragment in source, f"{path}: missing dispatch binding: {fragment}")


def forbid_fragments(path: Path, fragments: tuple[str, ...]) -> None:
    source = compact(path)
    for fragment in fragments:
        require(fragment not in source, f"{path}: stale incumbent dispatch: {fragment}")


def main() -> int:
    repo = Path(__file__).resolve().parents[2]
    avx = repo / "blst" / "avx"

    require_fragments(avx / "common" / "constexpr_utils.h", (
        '#include "../../common/defines.h"',
        "struct uint128_t",
        "struct int128_t",
        "constexpr BLST_FORCEINLINE uint128_t mul_u64_u64",
        "std::is_constant_evaluated()",
        "_umul128(a, b, &high)",
        "static_assert(sizeof(int128_t) == 16);",
    ))
    forbid_fragments(avx / "common" / "constexpr_utils.h", (
        "using int128_t = int64_t;",
        "__forceinline",
        "preprocess.h",
        "constexpr int64_t ceil_div",
        "constexpr int ceil_log2",
        "get_mask(",
        "mod_inverse(",
        "read_384bit",
        "parse_chunks",
    ))
    require(not (avx / "cpu" / "reduction" / "montgomery.hpp").exists(),
            "obsolete generic Montgomery owner was restored")
    require(not (avx / "rns" / "elementwise_reduce.h").exists(),
            "obsolete generic elementwise reducer was restored")
    require(not (avx / "rns" / "inversion.h").exists(),
            "orphaned VROOM test-only inversion chain was restored")

    bounds = avx / "rns" / "bounds.h"
    require_fragments(bounds, (
        'static_assert(left_shift < 63,',
        'constexpr int128_t multiplier =',
        'int64_t{1} << static_cast<int>(left_shift)',
    ))
    forbid_fragments(bounds, (
        'constexpr auto operator-() const',
        'static_cast<int64_t>(1) << left_shift',
    ))

    bounded = avx / "rns" / "bounded_element.h"
    require_fragments(bounded, (
        'static_assert(element_bits > 0 && element_bits <= 52,',
        'static constexpr int64_t max_wide_bound = int64_t{1} << (64 - 52 - 1);',
        'static_assert(BoundsType::lower >= 0,',
        'static_assert(BoundsType::upper < max_wide_bound,',
        'static_assert(BoundsType::lower >= OtherBounds::upper,',
        'return Zero(AVXVector<limbs>(0));',
    ))
    forbid_fragments(bounded, (
        '#include <iostream>',
        '#include <string>',
        'to_unsigned_array()',
        'bool check_bounds(',
        'std::array<int128_t, limbs> to_array()',
    ))

    ring = avx / "rns" / "ring_element.h"
    require_fragments(ring, (
        '// SPDX-License-Identifier: MIT',
        'using rns_bounds_type = RNSBounds;',
        'template<int64_t new_lower, int64_t new_upper>',
        'return ExpandedElement<decltype(result_m1), RNSBounds>',
        'Exact-representation predicates only.',
        'overload is intentionally ordered:',
    ))
    forbid_fragments(ring, (
        'template<int new_lower, int new_upper>',
        'ExpandedElement<decltype(result_m1), decltype(rns_bounds)>',
        'from_constant(',
        'ElementWideZero(',
    ))

    require_fragments(avx / "rns" / "matrixnok_fp.h", (
        "using FpRing = FpRingCandidate;",
        "static constexpr int64_t min_tracked_rns_multiple = -1932;",
        "static constexpr int64_t max_tracked_rns_multiple = 2377;",
        "max_tracked_rns_multiple - min_tracked_rns_multiple;",
        "combined_rns_multiple_bits + 2 * modulus_bits <= 2 * rns_bits",
        "RNSBoundsType::lower >= min_tracked_rns_multiple",
        "typename Small::rns_bounds_type",
        "AVXVector<limbs>& out_high",
        "The three vector objects must be distinct.",
        "raw high/low radix-2^52 accumulators consumed by `from_vec`",
        "return multiply(left, right);",
        "return multiply_batch(left, right);",
    ))
    forbid_fragments(avx / "rns" / "matrixnok_fp.h", (
        '#include "../../common/bytes.h"',
        "from_hex(",
        "std::remove_cv_t<decltype(small.rns_bounds)>",
        "AVXVector<8>& out_high",
        "return batch_reduce_expand(left * right)[0];",
        "products[index] = ready<MAX_ADD>(left[index] * right[index]);",
    ))
    forbid_fragments(avx / "curve" / "scalar_mult.h", (
        'ring.from_hex(',
        'GLVMult(const G1Curve &curve, const Ring &ring)',
        'GLSMult(const Ring &ring)',
    ))
    require_fragments(avx / "fields" / "matrixnok_fp12.h", (
        "using Fp12Ring = Fp12RingCandidate;",
    ))
    require_fragments(avx / "curve" / "matrixnok_g1.h", (
        "using G1Curve = G1CurveCandidate;",
        "using G1Projective = G1ProjectiveCandidate;",
        "using G1Affine = G1AffineCandidate;",
        "using G1GLV = G1GLVCandidate;",
    ))
    require_fragments(avx / "curve" / "matrixnok_g2.h", (
        "using G2Fp2Ring = G2Fp2RingCandidate;",
        "using G2Curve = G2CurveCandidate;",
        "using G2Projective = G2ProjectiveCandidate;",
        "using G2Affine = G2AffineCandidate;",
        "using G2GLS = G2GLSCandidate;",
    ))
    require_fragments(avx / "pairing" / "matrixnok_pairing.h", (
        "using MillerEngine = MillerCandidate;",
        "using PairingInverter = PairingInverterCandidate;",
    ))

    point = avx / "curve" / "point.h"
    require_fragments(point, (
        "using CurveRing = matrixnok::FpRing;",
        "using CurveFp2Ring = matrixnok::G2Fp2Ring;",
        "using CurveG1 = matrixnok::G1Curve;",
        "using CurveG2 = matrixnok::G2Curve;",
        "using CurveG1GLV = matrixnok::G1GLV;",
        "using CurveG2GLS = matrixnok::G2GLS;",
        "inline constexpr CurveRing kCurveRing{};",
        "inline constexpr CurveFp2Ring kCurveFp2Ring{kCurveRing};",
        "inline constexpr CurveG1GLV kCurveG1Glv(",
        "inline constexpr CurveG2GLS kCurveG2Gls(",
        "static constexpr G1Point g{",
        "static constexpr G2Point g{",
    ))
    forbid_fragments(point, (
        "static const CurveRing ring;",
        "static const CurveFp2Ring fp2ring(",
        "static const CurveG1GLV mult(",
        "static const CurveG2GLS mult(",
        "scalar_fp_to_ring(blst::scalar::get_BLST_Aprime_E1())",
        "scalar_fp2_to_ring(blst::scalar::get_BLST_Aprime_E2())",
        "const blst::scalar::G1AffinePoint& scalar_g",
        "const blst::scalar::G2AffinePoint& scalar_g",
    ))

    gt = avx / "curve" / "GT.h"
    require_fragments(gt, (
        "using Ring = matrixnok::FpRing;",
        "using FP12 = matrixnok::Fp12Ring;",
        "using PairingFp2 = matrixnok::G2Fp2Ring;",
        "matrixnok::generated::gt_generator",
        "matrixnok::PairingInverter inverter;",
        "const matrixnok::MillerEngine miller(",
        "inline constexpr FP12 kFp12Ring{};",
        "inline constexpr GT kGtIdentity(",
    ))
    forbid_fragments(gt, (
        "static const Ring value;",
        "static const PairingFp2 value(",
        "static const FP12 value(",
    ))

    fft = avx / "fft" / "gt_context.h"
    require_fragments(fft, (
        "using Ring = matrixnok::FpRing;",
        "using Fp12Ring = matrixnok::Fp12Ring;",
        "inline static constexpr Ring kRing{};",
        "inline static constexpr Fp12Ring kFp12Ring{};",
    ))
    forbid_fragments(fft, (
        "static const Ring instance;",
        "static const Fp12Ring instance(",
    ))

    avx_c_api = repo / "blst" / "avx_c_api.cpp"
    require_fragments(avx_c_api, (
        "&AvxG1::generator()",
        "&AvxG2::generator()",
        "&blst::avx::kGtIdentity",
    ))
    forbid_fragments(avx_c_api, (
        "static const AvxG1 generator",
        "static const AvxG2 generator",
        "static const AvxGt one",
    ))

    for owner in (point, gt, fft):
        forbid_fragments(owner, (
            "BoundedRing<",
            "generated::fp52",
            "AvxGt52",
            "scalar_gt_vector_to_avx52",
            "avx52_vector_to_scalar_gt",
        ))

    print("validated MatrixNoK production AVX dispatch")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
