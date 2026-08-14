#!/usr/bin/env python3
# Copyright 2026 Category Labs
# SPDX-License-Identifier: Apache-2.0

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path


REQUIRED_CATEGORIES = {"success", "rejection", "boundary", "aliasing", "roundtrip"}


def load_rows(root: Path) -> list[dict[str, str]]:
    path = root / "fuzz/COVERAGE_PATHS.csv"
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    expected = ["target", "category", "name", "file", "marker"]
    if not rows or list(rows[0]) != expected:
        raise ValueError(f"{path}: invalid or empty coverage-path ledger")
    return rows


def check_config(root: Path, rows: list[dict[str, str]]) -> None:
    categories = {row["category"] for row in rows}
    if categories != REQUIRED_CATEGORIES:
        raise ValueError(
            "coverage categories differ: "
            f"expected {sorted(REQUIRED_CATEGORIES)}, found {sorted(categories)}"
        )
    seen: set[tuple[str, str]] = set()
    for row in rows:
        key = (row["target"], row["marker"])
        if key in seen:
            raise ValueError(f"duplicate coverage path: {key}")
        seen.add(key)
        source = root / row["file"]
        if not source.is_file():
            raise ValueError(f"coverage source does not exist: {row['file']}")
        marker = f"FUZZ_PATH: {row['marker']}"
        if source.read_text(encoding="utf-8").count(marker) != 1:
            raise ValueError(f"coverage marker must occur exactly once: {marker}")


def executed_count(output: str, marker: str) -> str | None:
    needle = f"FUZZ_PATH: {marker}"
    for line in output.splitlines():
        if needle not in line:
            continue
        match = re.match(r"\s*\d+\|\s*([^|]*)\|", line)
        if not match:
            return None
        count = match.group(1).strip().replace(",", "")
        if not count or count.startswith("0") or set(count) == {"#"}:
            return None
        return count
    return None


def check_target(
    root: Path,
    rows: list[dict[str, str]],
    target: str,
    llvm_cov: str,
    binary: Path,
    profile: Path,
) -> None:
    selected = [row for row in rows if row["target"] == target]
    if not selected:
        return
    by_file: dict[str, list[dict[str, str]]] = {}
    for row in selected:
        by_file.setdefault(row["file"], []).append(row)
    for relative, file_rows in by_file.items():
        source = root / relative
        command = [
            llvm_cov,
            "show",
            str(binary),
            f"-instr-profile={profile}",
            "-show-line-counts-or-regions",
            "--sources",
            str(source),
        ]
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        for row in file_rows:
            count = executed_count(result.stdout, row["marker"])
            if count is None:
                raise ValueError(
                    f"uncovered named path {target}:{row['category']}:{row['name']}"
                )
            print(f"covered {target}:{row['category']}:{row['name']} count={count}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--check-config", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--target")
    parser.add_argument("--llvm-cov")
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--profile", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    try:
        rows = load_rows(root)
        check_config(root, rows)
        if args.self_test:
            if executed_count("  7|      3|x // FUZZ_PATH: sample", "sample") != "3":
                raise ValueError("coverage parser rejected an executed path")
            if executed_count("  7|      0|x // FUZZ_PATH: sample", "sample") is not None:
                raise ValueError("coverage parser accepted an uncovered path")
        if args.target:
            if not all((args.llvm_cov, args.binary, args.profile)):
                raise ValueError("target checks require llvm-cov, binary, and profile")
            check_target(
                root, rows, args.target, args.llvm_cov,
                args.binary.resolve(), args.profile.resolve(),
            )
        elif not (args.check_config or args.self_test):
            raise ValueError("select --check-config, --self-test, or --target")
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"fuzz coverage-path check failed: {error}", file=sys.stderr)
        return 1
    print("fuzz coverage-path check passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
