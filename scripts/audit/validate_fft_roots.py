#!/usr/bin/env python3
"""Independently validate the compile-time BLS12-381 Fr FFT root table."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


MODULUS = int(
    "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001",
    16,
)
TWO_ADICITY = 32
MONTGOMERY_R = pow(2, 256, MODULUS)
ENTRY = re.compile(
    r"Fr\{\{\s*"
    r"(0x[0-9a-fA-F]+)ULL,\s*"
    r"(0x[0-9a-fA-F]+)ULL,\s*"
    r"(0x[0-9a-fA-F]+)ULL,\s*"
    r"(0x[0-9a-fA-F]+)ULL\s*"
    r"\}\},\s*//\s*2\^(\d+)"
)


def expected_roots() -> list[tuple[int, int, int, int]]:
    roots = [0] * (TWO_ADICITY + 1)
    roots[TWO_ADICITY] = pow(
        7,
        (MODULUS - 1) // (1 << TWO_ADICITY),
        MODULUS,
    )
    for index in range(TWO_ADICITY - 1, -1, -1):
        roots[index] = roots[index + 1] * roots[index + 1] % MODULUS

    encoded = [root * MONTGOMERY_R % MODULUS for root in roots]
    return [
        tuple((value >> (64 * limb)) & ((1 << 64) - 1) for limb in range(4))
        for value in encoded
    ]


def parse_table(source: str) -> list[tuple[int, tuple[int, int, int, int]]]:
    return [
        (int(match.group(5)), tuple(int(match.group(i), 16) for i in range(1, 5)))
        for match in ENTRY.finditer(source)
    ]


def validate(source: str) -> None:
    actual = parse_table(source)
    if len(actual) != TWO_ADICITY + 1:
        raise ValueError(f"expected 33 root entries, found {len(actual)}")

    expected = expected_roots()
    for position, (index, limbs) in enumerate(actual):
        if index != position:
            raise ValueError(f"entry {position} is labelled 2^{index}")
        if limbs != expected[position]:
            raise ValueError(
                f"2^{position} root mismatch: expected {expected[position]!r}, got {limbs!r}"
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "header",
        nargs="?",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "blst/fft/primitive_root.h",
    )
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    source = args.header.read_text(encoding="utf-8")
    validate(source)

    if args.self_test:
        corrupted = source.replace("0x00000001fffffffeULL", "0x00000001ffffffffULL", 1)
        try:
            validate(corrupted)
        except ValueError:
            pass
        else:
            raise AssertionError("validator accepted a deliberately corrupted root table")

    print(f"validated {TWO_ADICITY + 1} independently derived Fr FFT roots")


if __name__ == "__main__":
    main()
