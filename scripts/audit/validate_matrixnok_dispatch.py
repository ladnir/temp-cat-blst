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
        '#include "preprocess.h"',
        "struct uint128_t",
        "struct int128_t",
        "constexpr INLINE uint128_t mul_u64_u64",
        "std::is_constant_evaluated()",
        "_umul128(a, b, &high)",
        "static_assert(sizeof(int128_t) == 16);",
    ))
    forbid_fragments(avx / "common" / "constexpr_utils.h", (
        "using int128_t = int64_t;",
        "__forceinline",
    ))

    require_fragments(avx / "rns" / "matrixnok_fp.h", (
        "using FpRing = FpRingCandidate;",
        "return multiply(left, right);",
        "return multiply_batch(left, right);",
    ))
    forbid_fragments(avx / "rns" / "matrixnok_fp.h", (
        "return batch_reduce_expand(left * right)[0];",
        "products[index] = ready<MAX_ADD>(left[index] * right[index]);",
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
