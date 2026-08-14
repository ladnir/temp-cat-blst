#!/usr/bin/env python3
# Copyright 2026 Category Labs
# SPDX-License-Identifier: Apache-2.0

"""Enforce the fuzz adapter boundary without compiling production code."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path


SOURCE_SUFFIXES = {".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp"}
BLST_ALLOWED_INTERNAL = {
    "fuzz/internal/adapter.cpp",
    "fuzz/internal/adapter.h",
    "fuzz/internal/adapter_access.h",
}
INTERNAL_INCLUDE_ALLOWED = BLST_ALLOWED_INTERNAL | {
    "fuzz/internal/factory_access.h",
    "fuzz/internal/scalar_factory.h",
    "fuzz/tests/scalar_factory_tests.cpp",
    "fuzz/tests/adapter_factory_tests.cpp",
    "fuzz/targets/decode.cpp",
    "fuzz/targets/bytes.cpp",
}
IGNORED_PREFIXES = ("fuzz/inventory/",)
ALLOWED_FACTORY_ADAPTER_CALLS = {
    "as_candidate",
    "as_on_curve",
    "incomplete_add_pair",
    "make_g1_candidate",
    "make_g2_candidate",
    "make_pairing_state",
    "pairing_commit",
    "precompute_middle_product_b",
    "prepare_fft_domain",
    "reduce_fp2",
    "reduce_fp12",
    "refine_nonzero",
    "retained_non_subgroup_g1",
    "retained_non_subgroup_g2",
    "run_bytes_case",
    "run_decode_case",
}


@dataclass(frozen=True)
class Finding:
    path: str
    line: int
    message: str


def line_number(text: str, offset: int) -> int:
    return text.count("\n", 0, offset) + 1


def inspect_source(
    relative: str, text: str, known_adapter_calls: set[str] | None = None
) -> list[Finding]:
    if relative.startswith(IGNORED_PREFIXES):
        return []
    rules = [
        (
            re.compile(r"\breinterpret_cast\s*<"),
            "reinterpret_cast is forbidden in fuzz code",
        ),
        (
            re.compile(r"\b(?:memcpy|memmove)\s*\([^;\n]*DomainValue"),
            "raw-byte fabrication of a domain value is forbidden",
        ),
        (
            re.compile(r"\bbit_cast\s*<[^>]*DomainValue"),
            "bit_cast fabrication of a domain value is forbidden",
        ),
    ]
    if relative not in INTERNAL_INCLUDE_ALLOWED:
        rules.append(
            (
                re.compile(r'^\s*#\s*include\s*[<"]fuzz/internal/', re.MULTILINE),
                "internal fuzz headers are allowed only in reviewed adapter and factory tests",
            )
        )
    if relative not in BLST_ALLOWED_INTERNAL:
        rules.extend((
        (
            re.compile(r'^\s*#\s*include\s*[<"](?:bindings/)?blst(?:/|_|\.)', re.MULTILINE),
            "public BLST headers are allowed only in the internal adapter",
        ),
        (
            re.compile(r"(?<!adapter::)\bblst_[A-Za-z0-9_]+\s*\("),
            "direct C API calls are allowed only in the internal adapter",
        ),
        ))
        if relative not in {"fuzz/support/domain.h", "fuzz/support/storage.h"}:
            rules.append(
                (
                    re.compile(r"\bLiveMiddleProductHandle\s*[<{(]"),
                    "opaque live handles must come from the adapter factory",
                )
            )
    findings: list[Finding] = []
    for pattern, message in rules:
        for match in pattern.finditer(text):
            findings.append(Finding(relative, line_number(text, match.start()), message))
    if known_adapter_calls is not None:
        for match in re.finditer(r"\badapter::([A-Za-z_][A-Za-z0-9_]*)\s*\(", text):
            if match.group(1) not in known_adapter_calls:
                findings.append(
                    Finding(
                        relative,
                        line_number(text, match.start()),
                        "adapter call is absent from the fuzz API ledger and factory allowlist",
                    )
                )
    return findings


def repository_sources(root: Path) -> list[tuple[str, str]]:
    fuzz_root = root / "fuzz"
    result: list[tuple[str, str]] = []
    for path in sorted(fuzz_root.rglob("*")):
        if path.is_file() and path.suffix.lower() in SOURCE_SUFFIXES:
            relative = path.relative_to(root).as_posix()
            result.append((relative, path.read_text(encoding="utf-8")))
    return result


def self_test() -> list[str]:
    failures: list[str] = []
    cases = {
        "forbidden include": ('#include "bindings/blst_scalar.h"\n', True),
        "direct call": ("void f(){ blst_p1_double(nullptr, nullptr); }\n", True),
        "reinterpret cast": ("auto p = reinterpret_cast<int*>(x);\n", True),
        "ordinary driver": ('#include "fuzz/support/domain.h"\n', False),
    }
    for name, (source, should_fail) in cases.items():
        failed = bool(inspect_source("fuzz/example.cpp", source))
        if failed != should_fail:
            failures.append(f"self-test '{name}' produced failed={failed}")
    if inspect_source("fuzz/internal/adapter.cpp", "blst_p1_double(0, 0);"):
        failures.append("self-test rejected an internal adapter call")
    if not inspect_source(
        "fuzz/example.cpp", "adapter::blst_unknown();", {"blst_p1_double"}
    ):
        failures.append("self-test accepted an unclassified adapter call")
    if inspect_source(
        "fuzz/example.cpp", "adapter::blst_p1_double();", {"blst_p1_double"}
    ):
        failures.append("self-test rejected a classified adapter call")
    if inspect_source(
        "fuzz/example.cpp", "adapter::refine_nonzero();", {"refine_nonzero"}
    ):
        failures.append("self-test rejected a classified factory call")
    return failures


def ledger_c_symbols(root: Path) -> set[str]:
    path = root / "docs/audit/FUZZ_API_LEDGER.csv"
    with path.open("r", encoding="utf-8", newline="") as handle:
        return {
            row["symbol"]
            for row in csv.DictReader(handle)
            if row["surface"] == "c"
        }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()

    errors = self_test() if args.self_test else []
    known_adapter_calls = ledger_c_symbols(root) | ALLOWED_FACTORY_ADAPTER_CALLS
    findings: list[Finding] = []
    for relative, text in repository_sources(root):
        findings.extend(inspect_source(relative, text, known_adapter_calls))

    for error in errors:
        print(error, file=sys.stderr)
    for finding in findings:
        print(f"{finding.path}:{finding.line}: {finding.message}", file=sys.stderr)
    if errors or findings:
        return 1
    print(f"fuzz structural gate passed: {len(repository_sources(root))} source files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
