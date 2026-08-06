#!/usr/bin/env python3
# Copyright Supranational LLC and Category Labs
# Licensed under the Apache License, Version 2.0, see LICENSE for details.
# SPDX-License-Identifier: Apache-2.0

"""Symbolically verify the generated BLS12-381 Fp addition chains."""

from dataclasses import dataclass
from pathlib import Path
import re
import sys


@dataclass(frozen=True)
class Chain:
    path: str
    target: int
    squares: int
    multiplies: int
    input_slot: int
    provenance: str


CHAINS = (
    Chain(
        "blst/scalar/fields/recip-addchain.h",
        4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559785,
        379,
        81,
        1,
        "79400f963d1464c53995d388ef1dcc38edbb8ebf",
    ),
    Chain(
        "blst/scalar/fields/sqrt-addchain.h",
        1000602388805416848354447456433976039139220704984751971333014534031007912622709466110671907282253916009473568139946,
        377,
        80,
        13,
        "d02b0d86e25d700ce4f4cb0eccac1f743aea6986",
    ),
)

DEST = r"(t\[(\d+)\]|out)"
ASSIGN = re.compile(r"t\[(\d+)\]\s*=\s*inp\s*;")
SQUARE = re.compile(rf"Field::sqr\({DEST},\s*t\[(\d+)\]\)")
MULTIPLY = re.compile(rf"Field::mul\({DEST},\s*t\[(\d+)\],\s*t\[(\d+)\]\)")
SQUARES_MULTIPLY = re.compile(
    rf"Field::sqr_n_mul\({DEST},\s*t\[(\d+)\],\s*(\d+),\s*t\[(\d+)\]\)"
)
TRACE = re.compile(r"/\*\s*(\d+):\s*([0-9a-f]+)\s*\*/")


def publish(exponents: dict[int, int], destination: str, slot: str | None,
            value: int) -> int | None:
    if destination == "out":
        return value
    assert slot is not None
    exponents[int(slot)] = value
    return None


def verify(root: Path, chain: Chain) -> list[str]:
    path = root / chain.path
    text = path.read_text(encoding="utf-8")
    errors: list[str] = []
    exponents: dict[int, int] = {}
    final: int | None = None
    squares = 0
    multiplies = 0
    next_trace = 0

    structural_requirements = (
        "namespace blst::scalar",
        "namespace detail",
        "std::is_same_v<Field, Fp>",
        chain.provenance,
    )
    for required in structural_requirements:
        if required not in text:
            errors.append(f"{chain.path}: missing structural marker {required!r}")
    for forbidden in ("vec_copy", "RECIPROCAL_MOD_BLS12_381_P"):
        if forbidden in text:
            errors.append(f"{chain.path}: forbidden legacy text {forbidden!r}")

    for line_number, line in enumerate(text.splitlines(), 1):
        stripped = line.strip()
        if stripped.startswith("/*") or stripped.startswith("//"):
            continue

        value: int | None = None
        matched = False
        expected_trace: int | None = None

        match = ASSIGN.search(line)
        if match:
            slot = int(match.group(1))
            value = 1
            exponents[slot] = value
            matched = True
            expected_trace = next_trace
            if slot != chain.input_slot:
                errors.append(
                    f"{chain.path}:{line_number}: input copied to t[{slot}], "
                    f"expected t[{chain.input_slot}]"
                )

        match = SQUARE.search(line)
        if match:
            destination, slot, source = match.groups()
            source_index = int(source)
            if source_index not in exponents:
                errors.append(
                    f"{chain.path}:{line_number}: read-before-write t[{source_index}]"
                )
            else:
                value = exponents[source_index] * 2
                result = publish(exponents, destination, slot, value)
                if result is not None:
                    if final is not None:
                        errors.append(f"{chain.path}:{line_number}: multiple output writes")
                    final = result
            squares += 1
            matched = True
            expected_trace = next_trace

        match = MULTIPLY.search(line)
        if match:
            destination, slot, left, right = match.groups()
            inputs = (int(left), int(right))
            missing = [index for index in inputs if index not in exponents]
            for index in missing:
                errors.append(
                    f"{chain.path}:{line_number}: read-before-write t[{index}]"
                )
            if not missing:
                value = exponents[inputs[0]] + exponents[inputs[1]]
                result = publish(exponents, destination, slot, value)
                if result is not None:
                    if final is not None:
                        errors.append(f"{chain.path}:{line_number}: multiple output writes")
                    final = result
            multiplies += 1
            matched = True
            expected_trace = next_trace

        match = SQUARES_MULTIPLY.search(line)
        if match:
            destination, slot, base, count, multiplier = match.groups()
            base_index = int(base)
            multiplier_index = int(multiplier)
            count_value = int(count)
            inputs = (base_index, multiplier_index)
            missing = [index for index in inputs if index not in exponents]
            for index in missing:
                errors.append(
                    f"{chain.path}:{line_number}: read-before-write t[{index}]"
                )
            if not missing:
                value = (exponents[base_index] << count_value) + exponents[multiplier_index]
                result = publish(exponents, destination, slot, value)
                if result is not None:
                    if final is not None:
                        errors.append(f"{chain.path}:{line_number}: multiple output writes")
                    final = result
            squares += count_value
            multiplies += 1
            matched = True
            expected_trace = next_trace + count_value

        if "Field::" in line and not matched:
            errors.append(f"{chain.path}:{line_number}: unrecognized active field operation")

        if matched:
            trace = TRACE.search(line)
            if trace is None:
                errors.append(f"{chain.path}:{line_number}: operation lacks exponent trace")
            else:
                actual_index = int(trace.group(1))
                assert expected_trace is not None
                if actual_index != expected_trace:
                    errors.append(
                        f"{chain.path}:{line_number}: trace index {actual_index}, "
                        f"expected {expected_trace}"
                    )
                if value is not None and int(trace.group(2), 16) != value:
                    errors.append(
                        f"{chain.path}:{line_number}: trace exponent does not match operation"
                    )
                next_trace = expected_trace + 1

    if final != chain.target:
        errors.append(f"{chain.path}: final exponent does not match the pinned target")
    if squares != chain.squares:
        errors.append(
            f"{chain.path}: {squares} squarings, expected {chain.squares}"
        )
    if multiplies != chain.multiplies:
        errors.append(
            f"{chain.path}: {multiplies} multiplications, expected {chain.multiplies}"
        )
    if set(exponents) != set(range(16)):
        errors.append(f"{chain.path}: scratch-slot coverage changed: {sorted(exponents)}")

    if not errors:
        print(
            f"{chain.path}: verified target, {squares} squarings, "
            f"{multiplies} multiplications, {next_trace} traced operations"
        )
    return errors


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    errors = [error for chain in CHAINS for error in verify(root, chain)]
    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        return 1
    print("BLS12-381 Fp addition-chain verification passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
