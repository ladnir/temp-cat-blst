#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Category Labs
"""Independently regenerate and validate the imported MatrixNoK constants.

This script deliberately uses only the Python standard library and parses the
checked-in C++ header. It does not import VROOM's generator or duplicate any
generated table, so it can detect a stale, truncated, or mistranscribed table.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Sequence


FP12_STANDARD_VALUES = {
    "one_third": int("11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d555526a9ffffffffc71d", 16),
    "neg_two_thirds": int("11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d555526a9ffffffffc71c", 16),
    "frob1_2": int("1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac", 16),
    "frob1_3": int("1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaad", 16),
    "frob1_4": int("5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe", 16),
    "frob1_5": int("1904d3bf02bb0667c231beb4202c0d1f0fd603fd3cbd5f4f7b2443d784bab9c4f67ea53d63e7813d8d0775ed92235fb8", 16),
    "frob1_6": int("fc3e2b36c4e03288e9e902231f9fb854a14787b6c7b36fec0c8ec971f63c5f282d5ac14d6c7ec22cf78a126ddc4af3", 16),
    "frob1_7": int("6af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09", 16),
    "frob1_8": int("135203e60180a68ee2e9c448d77a2cd91c3dedd930b1cf60ef396489f61eb45e304466cf3e67fa0af1ee7b04121bdea2", 16),
    "frob1_9": int("5b2cfd9013a5fd8df47fa6b48b1e045f39816240c0b8fee8beadf4d8e9c0566c63a3e6e257f87329b18fae980078116", 16),
    "frob1_10": int("144e4211384586c16bd3ad4afa99cc9170df3560e77982d0db45f3536814f0bd5871c1908bd478cd1ee605167ff82995", 16),
    "frob2_3": int("5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffeffff", 16),
}
FP12_STANDARD_ARRAY_NAMES = tuple(
    f"fp12_{name}_{basis}"
    for name in FP12_STANDARD_VALUES
    for basis in ("m1", "m2")
)


ARRAY_NAMES = (
    "moduli1",
    "moduli2",
    "moduli1_padded",
    "moduli2_padded",
    "one_m1",
    "one_m2",
    "modulus_in_rns_m1",
    "modulus_in_rns_m2",
    "six_modulus_in_rns_m1",
    "six_modulus_in_rns_m2",
    "max_neg_offset_m1",
    "max_neg_offset_m2",
    "fp12_offset_m1_high",
    "fp12_offset_m1_low",
    "fp12_offset_m2_high",
    "fp12_offset_m2_low",
    "r1_perm_mat_fixed",
    "r1_rns_mat",
    "r1_correction",
    "r1_shifted_quotient_estimations",
    "r2_rns_mat",
    "r2_correction",
    "r2_shifted_quotient_estimations",
    "convert_to_rns_mat",
    "convert_to_correction",
    "convert_to_shifted_quotient_estimations",
    "convert_from_rns_mat",
    "convert_from_correction",
    "convert_from_shifted_quotient_estimations",
    "reducer1_mont_factor",
    "reducer2_mont_factor",
    "gt_generator",
) + FP12_STANDARD_ARRAY_NAMES

INTEGER_NAMES = (
    "BLS12_381_RNS_DIM",
    "BLS12_381_PADDED_DIM",
    "BLS12_381_CONVERT_TO_DIGITS",
    "BLS12_381_CONVERT_TO_RADIX_BITS",
)

BIG_INTEGER_NAMES = (
    "target_modulus",
    "M1",
    "M2",
    "M1M2",
    "total_rotation",
    "total_inv_rotation",
    "qr_premult",
    "inv_factor_target",
    "r1_premult",
    "r1_postmult",
    "r2_premult",
    "r2_postmult",
    "convert_to_premult",
    "convert_to_postmult",
    "convert_from_premult",
    "convert_from_postmult",
)

GT_GENERATOR_SERIALIZED_HEX = (
    "1250ebd871fc0a92a7b2d83168d0d727272d441befa15c503dd8e90ce98db3e7"
    "b6d194f60839c508a84305aaca1789b6089a1c5b46e5110b86750ec6a53234886"
    "8a84045483c92b7af5af689452eafabf1a8943e50439f1d59882a98eaa0170f19f"
    "26337d205fb469cd6bd15c3d5a04dc88784fbb3d0b2dbdea54d43b2b73f2cbb12"
    "d58386a8703e0f948226e47ee89d06fba23eb7c5af0d9f80940ca771b6ffd5857b"
    "aaf222eb95a7d2809d61bfe02e1bfd1b68ff02f0b8102ae1c2d5d5ab1a1368bb4"
    "45c7c2d209703f239689ce34c0378a68e72a6b3b216da0e22a5031b54ddff57309"
    "396b38c881c4c849ec23e87193502b86edb8857c273fa075a50512937e0794e1e6"
    "5a7617c90d8bd66065b1fffe51d7a579973b1315021ec3c19934f11b8b424cd48b"
    "f38fcef68083b0b0ec5c81a93b330ee1a677d0d15ff7b984e8978ef48881e32fa"
    "c91b93b47333e2ba5703350f55a7aefcd3c31b4fcb6ce5771cc6a0e9786ab5973"
    "320c806ad360829107ba810c5a09ffdd9be2291a0c25a99a201b2f522473d17139"
    "1125ba84dc4007cfbf2f8da752f7c74185203fcca589ac719c34dffbbaad8431da"
    "d1c1fb597aaa5018107154f25a764bd3c79937a45b84546da634b8f6be14a8061"
    "e55cceba478b23f7dacaa35c8ca78beae9624045b4b604c581234d086a9902249b"
    "64728ffd21a189e87935a954051c7cdba7b3872629a4fafc05066245cb9108f024"
    "2d0fe3ef0f41e58663bf08cf068672cbd01a7ec73baca4d72ca93544deff686bfd"
    "6df543d48eaa24afe47e1efde449383b676631"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def parse_header(path: Path) -> tuple[dict[str, int], dict[str, list[int]]]:
    text = path.read_text(encoding="ascii")
    scalars: dict[str, int] = {}
    arrays: dict[str, list[int]] = {}

    for name in (
        "moduli1",
        "moduli2",
        "one_m1",
        "one_m2",
        "modulus_in_rns_m1",
        "modulus_in_rns_m2",
        "six_modulus_in_rns_m1",
        "six_modulus_in_rns_m2",
        "max_neg_offset_m1",
        "max_neg_offset_m2",
        "fp12_offset_m1_high",
        "fp12_offset_m1_low",
        "fp12_offset_m2_high",
        "fp12_offset_m2_low",
        "r1_perm_mat_fixed",
        "r2_rns_mat",
        "r2_correction",
        "convert_to_rns_mat",
        "convert_to_correction",
        "convert_from_rns_mat",
        "convert_from_correction",
        "gt_generator",
    ) + FP12_STANDARD_ARRAY_NAMES:
        require(
            re.search(rf"alignas\s*\(\s*64\s*\).*?\b{name}\s*=", text) is not None,
            f"{name} must retain 64-byte AVX row alignment",
        )

    for name in INTEGER_NAMES:
        match = re.search(rf"\b{name}\s*=\s*(\d+)\s*;", text)
        require(match is not None, f"missing integer {name}")
        scalars[name] = int(match.group(1))

    for name in BIG_INTEGER_NAMES:
        match = re.search(rf'\b{name}\s*\[\]\s*=\s*"(\d+)"\s*;', text)
        require(match is not None, f"missing big integer {name}")
        scalars[name] = int(match.group(1))

    for name in ARRAY_NAMES:
        match = re.search(rf"\b{name}\s*=\s*(.*?)\s*;", text, re.DOTALL)
        require(match is not None, f"missing array {name}")
        arrays[name] = [int(value) for value in re.findall(r"\d+", match.group(1))]

    return scalars, arrays


def product(values: Sequence[int]) -> int:
    result = 1
    for value in values:
        result *= value
    return result


def ceil_div(numerator: int, denominator: int) -> int:
    return numerator // denominator + (numerator % denominator != 0)


def words(value: int, bits: int, count: int) -> list[int]:
    require(0 <= value < (1 << (bits * count)), "radix decomposition overflow")
    mask = (1 << bits) - 1
    return [(value >> (bits * index)) & mask for index in range(count)]


def icrt_factors(moduli: Sequence[int]) -> list[int]:
    modulus = product(moduli)
    result = []
    for lane_modulus in moduli:
        rest = modulus // lane_modulus
        result.append(rest * pow(rest % lane_modulus, -1, lane_modulus) % modulus)
    return result


def regenerate_matrix(
    moduli_in: Sequence[int],
    moduli_out: Sequence[int],
    premult: int,
    postmult: int,
    *,
    input_digits: int,
    input_bits: int,
    output_digits: int,
    output_bits: int,
    round_quotient_down: bool = False,
) -> tuple[list[int], list[int], list[int]]:
    """Regenerate a mathematical wide-conversion table independently.

    VROOM rounds every fixed-point quotient coefficient upward. CatBLST keeps
    that policy for R1, R2, and ConvertFrom, but rounds ConvertTo downward so a
    canonical integer import cannot select a negative CRT representative.
    """
    modulus_in = product(moduli_in)
    modulus_out = product(moduli_out)
    fixed_point = 1 << 64
    factors = icrt_factors(moduli_in)
    matrix: list[int] = []
    estimates: list[int] = []

    # The exported table is input-modulus-major, then input-digit-major.
    for input_index in range(len(moduli_in)):
        for digit in range(input_digits):
            shifted = premult * (1 << (input_bits * digit)) * factors[input_index]
            shifted %= modulus_in
            target = shifted * postmult % modulus_out
            row: list[int] = []
            for lane_modulus in moduli_out:
                row.extend(words(target % lane_modulus, output_bits, output_digits))
            matrix.extend(row)
            scaled = shifted * fixed_point
            estimates.append(
                scaled // modulus_in
                if round_quotient_down
                else ceil_div(scaled, modulus_in)
            )

    correction_value = ((-modulus_in) % modulus_out) * postmult % modulus_out
    correction: list[int] = []
    for lane_modulus in moduli_out:
        correction.extend(words(correction_value % lane_modulus, output_bits, output_digits))
    return matrix, correction, estimates


def validate_shape(name: str, values: Sequence[int], expected: int) -> None:
    require(len(values) == expected, f"{name}: expected {expected} entries, found {len(values)}")


def validate_constants(s: dict[str, int], a: dict[str, list[int]]) -> None:
    require(s["BLS12_381_RNS_DIM"] == 8, "RNS dimension is not 8")
    require(s["BLS12_381_PADDED_DIM"] == 8, "padded RNS dimension is not 8")
    require(s["BLS12_381_CONVERT_TO_DIGITS"] == 8, "integer input does not have 8 digits")
    require(s["BLS12_381_CONVERT_TO_RADIX_BITS"] == 50, "integer input radix is not 50 bits")

    m1_moduli = a["moduli1"]
    m2_moduli = a["moduli2"]
    validate_shape("moduli1", m1_moduli, 8)
    validate_shape("moduli2", m2_moduli, 8)
    require(a["moduli1_padded"] == m1_moduli, "moduli1 padding changed productive lanes")
    require(a["moduli2_padded"] == m2_moduli, "moduli2 padding changed productive lanes")

    all_moduli = m1_moduli + m2_moduli
    require(len(set(all_moduli)) == 16, "RNS moduli are not distinct")
    for modulus in all_moduli:
        require(modulus & 1 == 1, f"RNS modulus {modulus} is even")
        require(modulus < 1 << 50, f"RNS modulus {modulus} exceeds 50 bits")
    for left_index, left in enumerate(all_moduli):
        for right in all_moduli[left_index + 1 :]:
            require(math.gcd(left, right) == 1, f"RNS moduli {left} and {right} are not coprime")

    p = s["target_modulus"]
    m1 = s["M1"]
    m2 = s["M2"]
    combined = s["M1M2"]
    require(m1 == product(m1_moduli), "M1 is not the product of moduli1")
    require(m2 == product(m2_moduli), "M2 is not the product of moduli2")
    require(combined == m1 * m2, "M1M2 is not M1 * M2")
    require(m1 > 6 * p and m2 > 6 * p, "RNS dynamic range is below the imported design bound")
    require(s["total_rotation"] * s["total_inv_rotation"] % combined == 1, "rotation inverse is invalid")
    require(s["inv_factor_target"] == pow(m1, -1, p), "target inverse factor is invalid")

    # VROOM's BLS12-381 schedule tracks p-multiples in [-1932, 2377].
    # Completion adds 1932*p, so the complete interval is [0, 4309].
    # Its wide redundant encoding must fit strictly inside M1*M2.
    min_tracked_rns = -1932
    max_tracked_rns = 2377
    completion_offset = -min_tracked_rns
    max_completed_rns = max_tracked_rns - min_tracked_rns
    base_redundant_multiple = 40
    require(completion_offset == 1932, "negative RNS completion offset changed")
    require(max_completed_rns == 4309, "completed RNS range changed")
    require(
        max_completed_rns * base_redundant_multiple**2 * p**2 < combined,
        "completed wide RNS range exceeds M1*M2",
    )

    def exact_montgomery(value: int, *, wide: bool = False) -> tuple[list[int], list[int]]:
        redundance = value // p
        encoded = value * m1 % p
        if wide:
            encoded = encoded * m1 % p
        encoded += redundance * p
        encoded = encoded * s["total_rotation"] % combined
        if wide:
            encoded = encoded * s["total_rotation"] % combined
        return (
            [encoded % modulus for modulus in m1_moduli],
            [encoded % modulus for modulus in m2_moduli],
        )

    fixed_values = (
        ("one", 1, False),
        ("modulus_in_rns", 40 * p, False),
        ("six_modulus_in_rns", 6 * 40 * p, False),
        ("max_neg_offset", completion_offset * base_redundant_multiple**2 * p * p, True),
    )
    for name, value, wide in fixed_values:
        expected_m1, expected_m2 = exact_montgomery(value, wide=wide)
        require(a[f"{name}_m1"] == expected_m1, f"{name} M1 encoding is invalid")
        require(a[f"{name}_m2"] == expected_m2, f"{name} M2 encoding is invalid")

    serialized_coefficients = [
        int(GT_GENERATOR_SERIALIZED_HEX[offset : offset + 96], 16)
        for offset in range(0, len(GT_GENERATOR_SERIALIZED_HEX), 96)
    ]
    require(len(serialized_coefficients) == 12,
        "GT generator serialization must contain twelve Fp coefficients")
    internal_order = (0, 1, 4, 5, 8, 9, 2, 3, 6, 7, 10, 11)
    expected_gt: list[int] = []
    for serialized_index in internal_order:
        expected_m1, expected_m2 = exact_montgomery(
            serialized_coefficients[serialized_index])
        expected_gt.extend(expected_m1)
        expected_gt.extend(expected_m2)
    require(a["gt_generator"] == expected_gt,
        "GT generator MatrixNoK encoding is invalid")

    for name, value in FP12_STANDARD_VALUES.items():
        expected_m1, expected_m2 = exact_montgomery(value)
        require(
            a[f"fp12_{name}_m1"] == expected_m1,
            f"Fp12 {name} M1 encoding is invalid",
        )
        require(
            a[f"fp12_{name}_m2"] == expected_m2,
            f"Fp12 {name} M2 encoding is invalid",
        )

    fp12_element_bound = 80
    fp12_rns_bound = 1440
    desired_m1, desired_m2 = exact_montgomery(
        9 * fp12_rns_bound * p * p, wide=True)
    fp12_minimum_low = (1 << 52) * fp12_element_bound
    fp12_mask = (1 << 52) - 1

    for base_name, moduli, desired in (
        ("m1", m1_moduli, desired_m1),
        ("m2", m2_moduli, desired_m2),
    ):
        expected_high: list[int] = []
        expected_low: list[int] = []
        for modulus, desired_residue in zip(moduli, desired):
            wide_product = fp12_element_bound * modulus * modulus
            high = wide_product >> 52
            low = (wide_product & fp12_mask) + desired_residue
            low += modulus * fp12_element_bound
            while low < fp12_minimum_low:
                low += modulus
            expected_high.append(high)
            expected_low.append(low)
        require(
            a[f"fp12_offset_{base_name}_high"] == expected_high,
            f"Fp12 {base_name.upper()} high offset is invalid",
        )
        require(
            a[f"fp12_offset_{base_name}_low"] == expected_low,
            f"Fp12 {base_name.upper()} low offset is invalid",
        )

    radix = 1 << 52
    inverse_factor = pow(m1, -1, m2)
    original_r1_factor = pow(-p, -1, m1) * pow(radix, -1, m1) % m1
    require(
        s["qr_premult"] * s["qr_premult"] * s["r1_premult"] % m1 == original_r1_factor,
        "quadratic-residue premultiplier relation is invalid",
    )
    require(s["r1_postmult"] == p * inverse_factor * inverse_factor * radix * radix, "r1 postmult is invalid")
    require(s["r2_premult"] == m1 * pow(radix, -1, m2), "r2 premult is invalid")
    require(s["r2_postmult"] == radix * radix * s["qr_premult"], "r2 postmult is invalid")
    require(s["convert_to_premult"] == m1, "convert-to premult is invalid")
    require(s["convert_to_postmult"] == inverse_factor * radix * radix, "convert-to postmult is invalid")
    require(s["convert_from_premult"] == s["r2_premult"], "convert-from premult is invalid")
    require(s["convert_from_postmult"] == pow(m1, -1, p), "convert-from postmult is invalid")

    icrt1 = m2 * pow(m2, -1, m1)
    icrt2 = m1 * pow(m1, -1, m2)
    rotation_factor = (s["qr_premult"] * icrt1 + inverse_factor * icrt2) % combined
    require(s["total_rotation"] == radix * rotation_factor % combined, "total rotation is invalid")

    for modulus, factor in zip(m1_moduli, a["reducer1_mont_factor"]):
        require(modulus * factor % radix == 1, f"invalid M1 Montgomery factor for {modulus}")
    for modulus, factor in zip(m2_moduli, a["reducer2_mont_factor"]):
        require(modulus * factor % radix == 1, f"invalid M2 Montgomery factor for {modulus}")

    cases = (
        ("r1", m1_moduli, m2_moduli, s["r1_premult"], s["r1_postmult"], 1, 52, 1, 52, False),
        ("r2", m2_moduli, m1_moduli, s["r2_premult"], s["r2_postmult"], 1, 52, 1, 52, False),
        (
            "convert_to",
            [p],
            m2_moduli,
            s["convert_to_premult"],
            s["convert_to_postmult"],
            8,
            50,
            1,
            52,
            True,
        ),
        (
            "convert_from",
            m2_moduli,
            [p],
            s["convert_from_premult"],
            s["convert_from_postmult"],
            1,
            52,
            8,
            52,
            False,
        ),
    )
    for (
        name,
        moduli_in,
        moduli_out,
        premult,
        postmult,
        in_digits,
        in_bits,
        out_digits,
        out_bits,
        round_quotient_down,
    ) in cases:
        matrix, correction, estimates = regenerate_matrix(
            moduli_in,
            moduli_out,
            premult,
            postmult,
            input_digits=in_digits,
            input_bits=in_bits,
            output_digits=out_digits,
            output_bits=out_bits,
            round_quotient_down=round_quotient_down,
        )
        require(a[f"{name}_rns_mat"] == matrix, f"{name} matrix does not regenerate")
        require(a[f"{name}_correction"] == correction, f"{name} correction does not regenerate")
        require(
            a[f"{name}_shifted_quotient_estimations"] == estimates,
            f"{name} quotient estimates do not regenerate",
        )

    base = a["r1_rns_mat"]
    permuted = a["r1_perm_mat_fixed"]
    validate_shape("r1_rns_mat", base, 64)
    validate_shape("r1_perm_mat_fixed", permuted, 64)
    expected_permuted = [base[((column - diagonal) % 8) * 8 + column] for diagonal in range(8) for column in range(8)]
    require(permuted == expected_permuted, "r1 cyclic permutation table does not match its base matrix")


def validate_rejects_corruption(scalars: dict[str, int], arrays: dict[str, list[int]]) -> None:
    """Demonstrate that representative single-entry corruptions fail closed."""
    mutations = (
        ("product", "scalar", "M1", 0),
        ("rotation", "scalar", "total_rotation", 0),
        ("matrix", "array", "r1_rns_mat", 17),
        ("permutation", "array", "r1_perm_mat_fixed", 31),
        ("quotient estimate", "array", "convert_to_shifted_quotient_estimations", 3),
        ("Montgomery factor", "array", "reducer2_mont_factor", 6),
        ("GT generator", "array", "gt_generator", 73),
    )
    for label, kind, name, index in mutations:
        changed_scalars = dict(scalars)
        changed_arrays = {key: list(value) for key, value in arrays.items()}
        if kind == "scalar":
            changed_scalars[name] ^= 1
        else:
            changed_arrays[name][index] ^= 1
        try:
            validate_constants(changed_scalars, changed_arrays)
        except ValueError:
            continue
        raise ValueError(f"validator accepted corrupted {label}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    default_header = Path(__file__).resolve().parents[2] / "blst" / "avx" / "generated_matrixnok_constants.hpp"
    parser.add_argument("header", nargs="?", type=Path, default=default_header)
    parser.add_argument("--self-test", action="store_true", help="also verify representative corruptions are rejected")
    args = parser.parse_args()
    scalars, arrays = parse_header(args.header)
    validate_constants(scalars, arrays)
    if args.self_test:
        validate_rejects_corruption(scalars, arrays)
    print(f"validated MatrixNoK constants: {args.header}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
