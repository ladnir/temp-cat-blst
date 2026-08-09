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

    require_fragments(avx / "rns" / "matrixnok_fp.h", (
        "using FpRing = FpRingCandidate;",
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
    ))

    gt = avx / "curve" / "GT.h"
    require_fragments(gt, (
        "using Ring = matrixnok::FpRing;",
        "using FP12 = matrixnok::Fp12Ring;",
        "using PairingFp2 = matrixnok::G2Fp2Ring;",
        "matrixnok::generated::gt_generator",
        "matrixnok::PairingInverter inverter;",
        "const matrixnok::MillerEngine miller(",
    ))

    fft = avx / "fft" / "gt_context.h"
    require_fragments(fft, (
        "using Ring = matrixnok::FpRing;",
        "using Fp12Ring = matrixnok::Fp12Ring;",
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
