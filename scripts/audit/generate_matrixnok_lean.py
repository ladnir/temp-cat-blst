#!/usr/bin/env python3
# Copyright 2026 Category Labs, Inc.
# SPDX-License-Identifier: Apache-2.0

"""Generate the Lean snapshot of the checked-in MatrixNoK C++ constants.

The source header is validated independently before emission. This script only
binds the exact reviewed C++ values into Lean; it does not prove their
mathematical properties. Those proofs live in `proof/avx/CatBLST/AVX/MatrixNoK`.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import re

import validate_matrixnok_constants as validator


MATRIX_NAMES = {
    "r1_perm_mat_fixed",
    "r1_rns_mat",
    "r2_rns_mat",
    "convert_to_rns_mat",
    "convert_from_rns_mat",
}


def lean_list(values: list[int], indent: str = "  ") -> str:
    return "[\n" + "\n".join(
        f"{indent}{value}{',' if index + 1 < len(values) else ''}"
        for index, value in enumerate(values)
    ) + "\n]"


def lean_matrix(values: list[int]) -> str:
    if len(values) != 64:
        raise ValueError(f"expected an 8x8 matrix, found {len(values)} entries")
    rows = [values[offset : offset + 8] for offset in range(0, 64, 8)]
    rendered_rows = []
    for row in rows:
        rendered_rows.append("  [" + ", ".join(str(value) for value in row) + "]")
    return "[\n" + ",\n".join(rendered_rows) + "\n]"


def render(header: Path) -> str:
    scalars, arrays = validator.parse_header(header)
    validator.validate_constants(scalars, arrays)

    text = header.read_text(encoding="ascii")
    skip_match = re.search(r"\br1_skip_correction\s*=\s*(true|false)\s*;", text)
    if skip_match is None:
        raise ValueError("missing r1_skip_correction")

    # Hash the parsed text, not the checkout's byte representation. Git may
    # materialize this header with either LF or CRLF line endings, but the
    # generated Lean constants must be identical on every host.
    source_hash = hashlib.sha256(text.encode("ascii")).hexdigest()
    lines = [
        "/-",
        "Copyright 2026 Category Labs, Inc.",
        "SPDX-License-Identifier: Apache-2.0",
        "",
        "Generated from blst/avx/generated_matrixnok_constants.hpp.",
        f"Source SHA-256: {source_hash}",
        "Do not edit by hand; run scripts/audit/generate_matrixnok_lean.py.",
        "-/",
        "",
        "namespace CatBLST.AVX.MatrixNoK.Generated",
        "",
        f'def sourceSha256 : String := "{source_hash}"',
        "",
    ]

    for name in validator.INTEGER_NAMES + validator.BIG_INTEGER_NAMES:
        lines.append(f"def {name} : Nat := {scalars[name]}")
    lines.extend([
        "",
        f"def r1_skip_correction : Bool := {skip_match.group(1)}",
        "",
    ])

    for name in validator.ARRAY_NAMES:
        values = arrays[name]
        if name in MATRIX_NAMES:
            lines.append(f"def {name} : List (List Nat) := {lean_matrix(values)}")
        else:
            lines.append(f"def {name} : List Nat := {lean_list(values)}")
        lines.append("")

    lines.append("end CatBLST.AVX.MatrixNoK.Generated")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--header",
        type=Path,
        default=repo / "blst" / "avx" / "generated_matrixnok_constants.hpp",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=repo / "proof" / "avx" / "CatBLST" / "AVX" / "MatrixNoK" / "Generated.lean",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail unless the checked-in Lean snapshot exactly matches",
    )
    args = parser.parse_args()

    expected = render(args.header.resolve())
    output = args.output.resolve()
    if args.check:
        if not output.exists() or output.read_text(encoding="ascii") != expected:
            raise ValueError(f"generated Lean constants are stale: {output}")
        print(f"validated generated Lean constants: {output}")
        return 0

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(expected, encoding="ascii", newline="\n")
    print(f"wrote generated Lean constants: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
