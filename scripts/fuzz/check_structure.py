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
    "fuzz/internal/adapter_avx.cpp",
    "fuzz/internal/adapter.h",
    "fuzz/internal/adapter_access.h",
    "fuzz/internal/adapter_common.h",
    "fuzz/internal/adapter_fields.cpp",
    "fuzz/internal/adapter_fft_middle.cpp",
    "fuzz/internal/adapter_groups.cpp",
    "fuzz/internal/adapter_protocol.cpp",
    "fuzz/internal/adapter_bls_api.cpp",
    "fuzz/internal/adapter_pairing.cpp",
    "fuzz/internal/adapter_msm.cpp",
    "fuzz/internal/adapter_uniq.cpp",
}
INTERNAL_INCLUDE_ALLOWED = BLST_ALLOWED_INTERNAL | {
    "fuzz/internal/factory_access.h",
    "fuzz/internal/scalar_factory.h",
    "fuzz/tests/scalar_factory_tests.cpp",
    "fuzz/tests/adapter_factory_tests.cpp",
    "fuzz/tests/seed_encoding_tests.cpp",
    "fuzz/targets/decode.cpp",
    "fuzz/targets/avx.cpp",
    "fuzz/targets/bytes.cpp",
    "fuzz/targets/groups.cpp",
    "fuzz/targets/fields.cpp",
    "fuzz/targets/fft_middle.cpp",
    "fuzz/targets/protocol.cpp",
    "fuzz/targets/bls_api.cpp",
    "fuzz/targets/pairing_state.cpp",
    "fuzz/targets/msm.cpp",
    "fuzz/targets/uniq.cpp",
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
    "run_avx_case",
    "run_abi_regressions",
    "run_decode_case",
    "run_fields_case",
    "run_fft_middle_case",
    "run_groups_case",
    "run_protocol_case",
    "run_bls_api_case",
    "run_pairing_state_case",
    "run_msm_case",
    "run_uniq_case",
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


def inspect_rust_source(relative: str, text: str) -> list[Finding]:
    findings: list[Finding] = []
    rules = (
        (re.compile(r"\bunsafe\s*\{"), "unsafe Rust is forbidden in safe fuzz targets"),
        (
            re.compile(r"\bblst_[a-z0-9_]+\s*\("),
            "raw Rust FFI calls are indirect through the C campaigns",
        ),
        (
            re.compile(r"\b(?:transmute|from_raw(?:_parts)?)\s*\("),
            "raw Rust object fabrication is forbidden in safe fuzz targets",
        ),
    )
    for pattern, message in rules:
        for match in pattern.finditer(text):
            findings.append(Finding(relative, line_number(text, match.start()), message))
    return findings


def inspect_rust_fuzz_targets(root: Path) -> list[Finding]:
    findings: list[Finding] = []
    target_root = root / "bindings/rust/fuzz/fuzz_targets"
    if not target_root.exists():
        return findings
    for path in sorted(target_root.glob("*.rs")):
        relative = path.relative_to(root).as_posix()
        findings.extend(
            inspect_rust_source(relative, path.read_text(encoding="utf-8"))
        )
    return findings


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
    if not inspect_rust_source("fuzz_targets/example.rs", "unsafe { raw(); }"):
        failures.append("self-test accepted unsafe Rust fuzz code")
    if not inspect_rust_source("fuzz_targets/example.rs", "blst_p1_double();"):
        failures.append("self-test accepted a raw Rust FFI call")
    if inspect_rust_source("fuzz_targets/example.rs", "blst_fr::default();"):
        failures.append("self-test rejected construction of a public Rust type")
    return failures


def ledger_c_symbols(root: Path) -> set[str]:
    path = root / "docs/audit/FUZZ_API_LEDGER.csv"
    with path.open("r", encoding="utf-8", newline="") as handle:
        return {
            row["symbol"]
            for row in csv.DictReader(handle)
            if row["surface"] == "c"
        }


def ledger_rows(root: Path) -> list[dict[str, str]]:
    path = root / "docs/audit/FUZZ_API_LEDGER.csv"
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def coverage_findings(root: Path) -> list[Finding]:
    rows = ledger_rows(root)
    findings: list[Finding] = []

    adapter_text = "\n".join(
        (root / relative).read_text(encoding="utf-8")
        for relative in sorted(BLST_ALLOWED_INTERNAL)
        if (root / relative).is_file()
    )
    required_c = {
        row["symbol"]
        for row in rows
        if row["surface"] == "c"
        and row["classification"] not in {"indirect", "excluded"}
    }
    for symbol in sorted(required_c):
        if not re.search(rf"\b{re.escape(symbol)}\s*\(", adapter_text):
            findings.append(
                Finding(
                    "docs/audit/FUZZ_API_LEDGER.csv",
                    1,
                    f"classified C API has no adapter call: {symbol}",
                )
            )

    required_cpp = {
        row["symbol"]
        for row in rows
        if row["surface"] == "cpp"
        and row["classification"] not in {"indirect", "excluded"}
    }
    claimed_cpp = set(
        re.findall(
            r"FUZZ_COVER:\s*([A-Za-z_][A-Za-z0-9_:~]+)",
            adapter_text,
        )
    )
    for symbol in sorted(required_cpp - claimed_cpp):
        findings.append(
            Finding(
                "docs/audit/FUZZ_API_LEDGER.csv",
                1,
                f"unique native C++ API has no exact coverage marker: {symbol}",
            )
        )
    for symbol in sorted(claimed_cpp - required_cpp):
        findings.append(
            Finding(
                "fuzz/internal",
                1,
                f"native C++ coverage marker is absent from the ledger: {symbol}",
            )
        )

    evidence_path = root / "bindings/rust/fuzz/SAFE_API_COVERAGE.txt"
    claimed_rust: dict[str, str] = {}
    if not evidence_path.is_file():
        findings.append(
            Finding(
                "bindings/rust/fuzz/SAFE_API_COVERAGE.txt",
                1,
                "safe Rust coverage ownership file is missing",
            )
        )
        return findings
    for number, raw in enumerate(
        evidence_path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.count("|") != 1:
            findings.append(
                Finding(evidence_path.relative_to(root).as_posix(), number,
                        "malformed safe Rust coverage row")
            )
            continue
        symbol, owner = line.split("|", 1)
        if symbol in claimed_rust:
            findings.append(
                Finding(evidence_path.relative_to(root).as_posix(), number,
                        f"duplicate safe Rust coverage row: {symbol}")
            )
            continue
        claimed_rust[symbol] = owner

    required_rust = {
        row["symbol"]
        for row in rows
        if row["surface"] == "rust_safe"
        and row["classification"] != "excluded"
    }
    for symbol in sorted(required_rust - claimed_rust.keys()):
        findings.append(
            Finding(evidence_path.relative_to(root).as_posix(), 1,
                    f"safe Rust API has no evidence owner: {symbol}")
        )
    for symbol in sorted(claimed_rust.keys() - required_rust):
        findings.append(
            Finding(evidence_path.relative_to(root).as_posix(), 1,
                    f"safe Rust evidence row is absent from the ledger: {symbol}")
        )
    fuzz_root = root / "bindings/rust/fuzz"
    for symbol, owner in sorted(claimed_rust.items()):
        owner_path = fuzz_root / owner
        if not owner_path.is_file():
            findings.append(
                Finding(evidence_path.relative_to(root).as_posix(), 1,
                        f"safe Rust evidence owner does not exist: {owner}")
            )
            continue
        leaf = symbol.rsplit("::", 1)[-1]
        owner_text = owner_path.read_text(encoding="utf-8")
        if not re.search(rf"\b{re.escape(leaf)}\s*\(", owner_text):
            findings.append(
                Finding(owner_path.relative_to(root).as_posix(), 1,
                        f"claimed safe Rust API spelling is not executed: {symbol}")
            )

    coverage_source = root / "bindings/rust/fuzz/fuzz_targets/api_coverage.rs"
    if coverage_source.is_file():
        coverage_text = coverage_source.read_text(encoding="utf-8")
        for variant in ("min_pk", "min_sig"):
            if f"exercise_variant!({variant}" not in coverage_text:
                findings.append(
                    Finding(coverage_source.relative_to(root).as_posix(), 1,
                            f"shared safe wrapper campaign omits {variant}")
                )
    return findings


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
    findings.extend(inspect_rust_fuzz_targets(root))
    findings.extend(coverage_findings(root))

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
