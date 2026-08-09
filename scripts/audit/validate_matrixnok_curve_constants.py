#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Category Labs
"""Independently regenerate the AVX MatrixNoK curve constants.

The curve header contains compile-time RNS encodings used to avoid guarded
first-use initialization in production. This validator starts from canonical
BLS12-381 constants and the independently validated MatrixNoK parameters; it
does not trust or execute the C++ conversion implementation.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


P = int(
    "1a0111ea397fe69a4b1ba7b6434bacd7"
    "64774b84f38512bf6730d2a0f6b0f624"
    "1eabfffeb153ffffb9feffffffffaaab",
    16,
)

CANONICAL_VALUES = {
    "g1_generator_x": int("17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb", 16),
    "g1_generator_y": int("08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1", 16),
    "g2_generator_x0": int("024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8", 16),
    "g2_generator_x1": int("13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e", 16),
    "g2_generator_y0": int("0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801", 16),
    "g2_generator_y1": int("0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be", 16),
    "aprime_e1": int("00144698a3b8e9433d693a02c96d4982b0ea985383ee66a8d8e8981aefd881ac98936f8da0e0f97f5cf428082d584c1d", 16),
    "aprime_e2_x": 0,
    "aprime_e2_y": 240,
}
CANONICAL_VALUES["g1_generator_neg_y"] = (-CANONICAL_VALUES["g1_generator_y"]) % P
CANONICAL_VALUES["g2_generator_neg_y0"] = (-CANONICAL_VALUES["g2_generator_y0"]) % P
CANONICAL_VALUES["g2_generator_neg_y1"] = (-CANONICAL_VALUES["g2_generator_y1"]) % P


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def parse_array(text: str, name: str) -> list[int]:
    require(
        re.search(rf"alignas\s*\(\s*64\s*\).*?\b{name}\s*=", text) is not None,
        f"{name} must retain 64-byte AVX row alignment",
    )
    match = re.search(rf"\b{name}\s*=\s*\{{([^}}]*)\}}\s*;", text)
    require(match is not None, f"missing array {name}")
    values = [int(value) for value in re.findall(r"\d+", match.group(1))]
    require(len(values) == 8, f"{name}: expected 8 lanes, found {len(values)}")
    return values


def parse_decimal(base_text: str, name: str) -> int:
    match = re.search(rf'\b{name}\s*\[\]\s*=\s*"(\d+)"\s*;', base_text)
    require(match is not None, f"missing MatrixNoK parameter {name}")
    return int(match.group(1))


def parse_base_array(base_text: str, name: str) -> list[int]:
    match = re.search(rf"\b{name}\s*=\s*\{{([^}}]*)\}}\s*;", base_text)
    require(match is not None, f"missing MatrixNoK array {name}")
    values = [int(value) for value in re.findall(r"\d+", match.group(1))]
    require(len(values) == 8, f"{name}: expected 8 lanes, found {len(values)}")
    return values


def expected_encodings(base_text: str, value: int) -> tuple[list[int], list[int]]:
    moduli1 = parse_base_array(base_text, "moduli1")
    moduli2 = parse_base_array(base_text, "moduli2")
    m1 = parse_decimal(base_text, "M1")
    combined = parse_decimal(base_text, "M1M2")
    rotation = parse_decimal(base_text, "total_rotation")
    require(parse_decimal(base_text, "target_modulus") == P, "base target modulus is not BLS12-381 Fp")

    encoded = (value % P) * m1 % P
    encoded = encoded * rotation % combined
    return (
        [encoded % modulus for modulus in moduli1],
        [encoded % modulus for modulus in moduli2],
    )


def validate_constants(base_text: str, curve_text: str) -> None:
    for name, value in CANONICAL_VALUES.items():
        expected_m1, expected_m2 = expected_encodings(base_text, value)
        actual_m1 = parse_array(curve_text, f"{name}_m1")
        actual_m2 = parse_array(curve_text, f"{name}_m2")
        require(actual_m1 == expected_m1, f"{name} M1 encoding is invalid")
        require(actual_m2 == expected_m2, f"{name} M2 encoding is invalid")


def validate_rejects_corruption(base_text: str, curve_text: str) -> None:
    match = re.search(r"(\bg1_generator_x_m1\s*=\s*\{\s*)(\d+)", curve_text)
    require(match is not None, "cannot construct curve-constant corruption self-test")
    corrupted = curve_text[: match.start(2)] + str(int(match.group(2)) ^ 1) + curve_text[match.end(2) :]
    try:
        validate_constants(base_text, corrupted)
    except ValueError:
        return
    raise ValueError("validator accepted a corrupted curve constant")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[2]
    parser.add_argument(
        "header",
        nargs="?",
        type=Path,
        default=root / "blst" / "avx" / "generated_matrixnok_curve_constants.hpp",
    )
    parser.add_argument(
        "--base-header",
        type=Path,
        default=root / "blst" / "avx" / "generated_matrixnok_constants.hpp",
    )
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    base_text = args.base_header.read_text(encoding="ascii")
    curve_text = args.header.read_text(encoding="ascii")
    validate_constants(base_text, curve_text)
    if args.self_test:
        validate_rejects_corruption(base_text, curve_text)
    print(f"validated MatrixNoK curve constants: {args.header}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
