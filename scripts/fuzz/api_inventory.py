#!/usr/bin/env python3
# Copyright 2026 Category Labs
# SPDX-License-Identifier: Apache-2.0

"""Generate and validate the symbol-level fuzz API ledger.

C declarations are parsed directly from their reviewed headers. Native C++
methods and safe Rust functions are consumed from Clang AST JSON and rustdoc
JSON, respectively. The generator does not infer mathematical validity from a
signature. It assigns a conservative campaign and domain description that a
reviewer must accept in the checked-in ledger.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator


LEDGER_COLUMNS = (
    "surface",
    "symbol",
    "owner",
    "campaign",
    "classification",
    "required_domain",
    "factory",
    "allowed_invalid_cases",
    "oracles",
    "notes",
)

ALLOWED_SURFACES = {"c", "cpp", "rust_safe", "rust_raw"}
ALLOWED_CAMPAIGNS = {
    "decode",
    "bytes",
    "fields",
    "groups",
    "protocol",
    "pairing_state",
    "uniq",
    "fft_middle",
    "msm",
    "avx",
    "regression",
    "excluded",
}
ALLOWED_CLASSIFICATIONS = {
    "direct",
    "constructed",
    "oracle",
    "indirect",
    "regression_only",
    "excluded",
}


@dataclass(frozen=True, order=True)
class ApiSymbol:
    surface: str
    symbol: str
    owner: str


def normalize_path(value: str, root: Path) -> str:
    path = Path(value.replace("\\", "/"))
    if path.is_absolute():
        try:
            path = path.resolve().relative_to(root.resolve())
        except ValueError:
            return path.as_posix()
    return path.as_posix().lstrip("./")


def strip_c_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    return re.sub(r"//[^\n]*", "", text)


def extract_c_symbols(root: Path, headers: Iterable[str]) -> list[ApiSymbol]:
    pattern = re.compile(
        r"(?ms)^\s*(?:BLST_[A-Z_]+\s+)*"
        r"[A-Za-z_][A-Za-z0-9_\s*]*?\b(blst_[A-Za-z0-9_]+)"
        r"\s*\([^;]*?\)\s*;"
    )
    result: list[ApiSymbol] = []
    for owner in headers:
        text = strip_c_comments((root / owner).read_text(encoding="utf-8"))
        names = sorted(set(pattern.findall(text)))
        result.extend(ApiSymbol("c", name, owner) for name in names)
    return result


def ast_location_file(node: dict[str, Any], inherited: str | None) -> str | None:
    loc = node.get("loc", {})
    direct = loc.get("file")
    if direct:
        return str(direct)
    begin = node.get("range", {}).get("begin", {})
    direct = begin.get("file")
    if direct:
        return str(direct)
    return inherited


def extract_cpp_symbols(
    root: Path,
    ast_path: Path,
    owners: Iterable[str],
    public_types: dict[str, str],
) -> list[ApiSymbol]:
    if ast_path.suffix == ".txt":
        return extract_cpp_symbols_from_text(root, ast_path, public_types)
    ast = json.loads(ast_path.read_text(encoding="utf-8"))
    owner_set = {normalize_path(owner, root) for owner in owners}
    result: set[ApiSymbol] = set()

    def visit(
        node: dict[str, Any],
        scope: tuple[str, ...],
        inherited_file: str | None,
        access: str = "public",
    ) -> None:
        kind = node.get("kind", "")
        name = node.get("name")
        current_file = ast_location_file(node, inherited_file)
        next_scope = scope

        if kind == "NamespaceDecl" and name:
            next_scope = scope + (str(name),)
        elif kind in {"CXXRecordDecl", "ClassTemplateDecl"} and name:
            next_scope = scope + (str(name),)

        if kind in {"FunctionDecl", "CXXMethodDecl", "CXXConstructorDecl"}:
            if name and access == "public" and current_file:
                owner = normalize_path(current_file, root)
                if owner in owner_set:
                    qualified = "::".join(scope + (str(name),))
                    if qualified.startswith("blst::"):
                        result.add(ApiSymbol("cpp", qualified, owner))

        child_access = access
        if kind == "CXXRecordDecl":
            child_access = "public" if node.get("tagUsed") == "struct" else "private"

        for child in node.get("inner", []):
            if child.get("kind") == "AccessSpecDecl":
                child_access = str(child.get("access", child_access))
                continue
            visit(child, next_scope, current_file, child_access)

    visit(ast, tuple(), None)
    return sorted(result)


def extract_cpp_symbols_from_text(
    root: Path, ast_path: Path, public_types: dict[str, str]
) -> list[ApiSymbol]:
    del root
    by_name: dict[str, list[tuple[str, str]]] = {}
    for qualified, owner in public_types.items():
        by_name.setdefault(qualified.split("::")[-1], []).append((qualified, owner))

    record_pattern = re.compile(
        r"CXXRecordDecl\s+0x[0-9a-f]+.*?\b(?:struct|class)\s+([A-Za-z_][A-Za-z0-9_]*)\s+definition"
    )
    callable_pattern = re.compile(
        r"(?:CXXMethodDecl|CXXConstructorDecl|CXXConversionDecl|CXXDestructorDecl)\s+"
        r"0x[0-9a-f]+.*?\b([A-Za-z_~][A-Za-z0-9_~]*|operator[^ '\"]*)\s+'"
    )

    stack: list[tuple[int, str, str]] = []
    result: set[ApiSymbol] = set()
    profile = "scalar"
    for line in ast_path.read_text(encoding="utf-8").splitlines():
        if line.startswith("# profile:"):
            profile = line.split(":", 1)[1].strip()
            stack.clear()
            continue
        kind_offsets = [offset for token in ("CXXRecordDecl", "CXXMethodDecl", "CXXConstructorDecl", "CXXConversionDecl", "CXXDestructorDecl") if (offset := line.find(token)) >= 0]
        if not kind_offsets:
            continue
        kind_offset = min(kind_offsets)
        # The textual tree uses two-character indentation cells. The byte
        # offset of the declaration kind therefore preserves parent/child
        # ordering even when the connector changes from `- to |-.
        depth = kind_offset
        while stack and stack[-1][0] >= depth:
            stack.pop()

        record_match = record_pattern.search(line)
        if record_match:
            candidates = by_name.get(record_match.group(1), [])
            if len(candidates) > 1 and profile in {"scalar", "avx"}:
                candidates = [
                    candidate
                    for candidate in candidates
                    if f"::{profile}::" in candidate[0]
                ]
            if len(candidates) == 1:
                qualified, owner = candidates[0]
                stack.append((depth, qualified, owner))
            else:
                # Preserve tree depth so methods of a nested implementation
                # type cannot be attributed to the nearest public ancestor.
                stack.append((depth, "", ""))
            continue

        callable_match = callable_pattern.search(line)
        if not callable_match or not stack:
            continue
        _, qualified, owner = stack[-1]
        if not qualified:
            continue
        name = callable_match.group(1)
        if " implicit " in f" {line} ":
            continue
        result.add(ApiSymbol("cpp", f"{qualified}::{name}", owner))
    return sorted(result)


def rustdoc_path(item_id: str, rustdoc: dict[str, Any]) -> str | None:
    paths = rustdoc.get("paths", {})
    path = paths.get(item_id)
    if not path:
        return None
    segments = path.get("path", [])
    return "::".join(str(segment) for segment in segments)


def rust_type_name(type_value: Any, rustdoc: dict[str, Any]) -> str:
    if not isinstance(type_value, dict):
        return "unknown"
    resolved = type_value.get("resolved_path")
    if isinstance(resolved, dict):
        item_id = resolved.get("id")
        path = rustdoc_path(str(item_id), rustdoc) if item_id is not None else None
        return path or str(resolved.get("name", "unknown"))
    generic = type_value.get("generic")
    if generic:
        return str(generic)
    borrowed = type_value.get("borrowed_ref")
    if isinstance(borrowed, dict):
        return rust_type_name(borrowed.get("type"), rustdoc)
    return "unknown"


def rust_visibility_is_public(item: dict[str, Any]) -> bool:
    visibility = item.get("visibility")
    return visibility == "public" or isinstance(visibility, dict) and "public" in visibility


def rust_function_is_safe(item: dict[str, Any]) -> bool:
    function = item.get("inner", {}).get("function")
    if not isinstance(function, dict):
        return False
    return not bool(function.get("header", {}).get("is_unsafe", False))


def extract_rust_symbols(
    root: Path,
    rustdoc_path_value: Path,
    manifest_owner: str,
    excluded_modules: Iterable[str],
) -> list[ApiSymbol]:
    del root
    rustdoc = json.loads(rustdoc_path_value.read_text(encoding="utf-8"))
    index: dict[str, dict[str, Any]] = rustdoc.get("index", {})
    excluded = tuple(excluded_modules)
    result: set[ApiSymbol] = set()

    def is_excluded(symbol: str) -> bool:
        return any(
            symbol == f"blst::{module}" or symbol.startswith(f"blst::{module}::")
            for module in excluded
        )

    for item_id, item in index.items():
        inner = item.get("inner", {})
        if (
            "function" in inner
            and rust_visibility_is_public(item)
            and rust_function_is_safe(item)
        ):
            symbol = rustdoc_path(item_id, rustdoc)
            if symbol and symbol.startswith("blst::") and not is_excluded(symbol):
                result.add(ApiSymbol("rust_safe", symbol, manifest_owner))

        impl = inner.get("impl")
        if not isinstance(impl, dict):
            continue
        target = rust_type_name(impl.get("for"), rustdoc)
        if not target.startswith("blst::") or is_excluded(target):
            continue
        for method_id in impl.get("items", []):
            method = index.get(str(method_id), {})
            if not rust_visibility_is_public(method):
                continue
            method_inner = method.get("inner", {})
            if "function" not in method_inner or not rust_function_is_safe(method):
                continue
            method_name = method.get("name")
            if method_name:
                result.add(
                    ApiSymbol(
                        "rust_safe",
                        f"{target}::{method_name}",
                        manifest_owner,
                    )
                )

    return sorted(result)


def campaign_for(api: ApiSymbol) -> str:
    name = api.symbol.lower()
    owner = api.owner.replace("\\", "/").lower()
    if api.surface == "rust_raw":
        return "avx" if "avx" in owner or "avx" in name else campaign_for(
            ApiSymbol("c", api.symbol.split("::")[-1], api.owner)
        )
    if "avx" in owner or "::avx" in name or "blst_avx_" in name:
        return "avx"
    if "uniq" in name:
        return "uniq"
    if any(token in owner for token in ("/fft/", "middle_product")) or any(
        token in name for token in ("fft", "interpol", "middle_product")
    ):
        return "fft_middle"
    if any(token in name for token in ("pairing", "miller", "finalverify", "precompute_lines")):
        return "pairing_state"
    if any(token in name for token in ("pippenger", "wbits", "mult_batch", "tile_", "mult_fr", "multpoint")):
        return "msm"
    if any(
        token in name
        for token in (
            "uncompress",
            "deserialize",
            "aggregate_in_g1",
            "aggregate_in_g2",
            "aggregated_in_g1",
            "aggregated_in_g2",
        )
    ):
        return "decode"
    if any(
        token in name
        for token in (
            "from_be",
            "from_le",
            "to_be",
            "to_le",
            "bendian",
            "lendian",
            "hexascii",
            "from_uint",
            "to_uint",
            "from_scalar",
            "to_scalar",
            "fr_check",
            "sk_check",
        )
    ):
        return "bytes"
    if any(
        token in name
        for token in (
            "key_gen",
            "keygen",
            "derive_",
            "sign",
            "verify",
            "hash_to",
            "encode_to",
            "map_to",
            "sha256",
            "expand_message",
            "sk_to_pk",
            "publickey",
            "secretkey",
            "signature",
            "aggregatepublickey",
            "aggregatesignature",
        )
    ):
        return "protocol"
    if any(
        token in name
        for token in ("fp12", "fp2", "fp", "fr", "gt", "field", "final_exp", "sk_add", "sk_sub", "sk_mul", "sk_inverse")
    ):
        return "fields"
    if any(token in name for token in ("p1", "p2", "g1", "g2", "point", "affine")):
        return "groups"
    return "regression"


def classification_for(api: ApiSymbol, campaign: str) -> str:
    name = api.symbol.lower()
    if api.surface == "rust_raw":
        return "indirect"
    if api.surface == "cpp":
        unique_methods = (
            "deserialize_unchecked",
            "dadd",
            "dadd_affine",
            "fft_forward_truncated_prefix",
            "fft_inverse_truncated_range",
        )
        return "constructed" if any(token in name for token in unique_methods) else "indirect"
    if any(token in name for token in ("scratch_sizeof", "scratch_len")):
        return "constructed"
    if any(
        token in name
        for token in (
            "uncompress",
            "deserialize",
            "aggregate_in_g1",
            "aggregate_in_g2",
        )
    ) and "aggregated_in_" not in name:
        return "direct"
    if any(token in name for token in ("sizeof", "alignof")):
        return "regression_only"
    if any(
        token in name
        for token in (
            "generator",
            "identity",
            "::one",
            "_one",
            "is_equal",
            "is_identity",
            "is_inf",
            "on_curve",
            "in_g1",
            "in_g2",
            "in_group",
            "_check",
            "domain_len",
            "scratch_sizeof",
            "get_dst",
        )
    ):
        return "oracle"
    if campaign == "decode":
        return "direct"
    if campaign == "bytes":
        if any(
            token in name
            for token in (
                "from_",
                "parse",
                "check",
                "hash_to",
            )
        ):
            return "direct"
        return "constructed"
    if campaign == "regression":
        return "regression_only"
    return "constructed"


DOMAIN_BY_CAMPAIGN = {
    "decode": "arbitrary bounded bytes and documented readable extent",
    "bytes": "arbitrary bounded bytes or words",
    "fields": "typed reduced field value; nonzero or GT refinement when required",
    "groups": "typed reduced coordinates, on-curve point, or subgroup point as declared",
    "protocol": "bounded byte slices or typed secret/subgroup operands",
    "pairing_state": "typed subgroup operands and valid pairing-state transition",
    "uniq": "owned bounded message slices and live uniqueness-tree storage",
    "fft_middle": "typed reduced vectors, prepared domain, or live middle-product handle",
    "msm": "bounded subgroup points, public scalars, and typed disjoint scratch",
    "avx": "typed scalar-equivalent value and supported AVX execution profile",
    "regression": "deterministic construction",
    "excluded": "outside supported fuzz scope",
}

FACTORY_BY_CAMPAIGN = {
    "decode": "bounded_bytes",
    "bytes": "bounded_bytes_or_words",
    "fields": "field_domain_factory",
    "groups": "group_domain_factory",
    "protocol": "protocol_input_factory",
    "pairing_state": "pairing_state_machine",
    "uniq": "uniq_call_frame",
    "fft_middle": "fft_middle_domain_factory",
    "msm": "msm_call_frame",
    "avx": "scalar_to_avx_domain_factory",
    "regression": "deterministic_fixture",
    "excluded": "not_applicable",
}

ORACLE_BY_CAMPAIGN = {
    "decode": "failure atomicity; curve predicate; canonical round trip",
    "bytes": "independent modulus oracle; endian symmetry; canonical round trip",
    "fields": "algebraic identities; canonical export; documented aliasing",
    "groups": "group identities; scalar oracle; serialization round trip",
    "protocol": "determinism; subgroup output; sign/verify agreement",
    "pairing_state": "state preservation; one-shot and incremental agreement",
    "uniq": "model set; capacity and duplicate result",
    "fft_middle": "round trip; bounded direct formula; failure atomicity",
    "msm": "direct sum of single-point multiplication; failure atomicity",
    "avx": "scalar differential; supported-profile sanitizer",
    "regression": "deterministic expected value or ABI property",
    "excluded": "not_applicable",
}


def required_domain_for(api: ApiSymbol, campaign: str, classification: str) -> str:
    name = api.symbol.lower()
    if classification == "indirect":
        return "inherits the domain of the corresponding direct owner"
    if campaign == "decode":
        return "arbitrary bounded bytes with the documented readable extent"
    if campaign == "bytes":
        if any(token in name for token in ("from_", "parse", "check", "hash_to")):
            return "arbitrary bounded bytes, words, or NUL-terminated hexadecimal text"
        return "typed scalar or reduced field input"
    if campaign == "fields":
        if "final_exp" in name:
            return "nonzero general Fp12 or Miller value"
        if "inverse" in name:
            return "Nonzero<ReducedField> of the declared field"
        if "cyclotomic" in name or "gt_" in name or "::gt::" in name:
            return "GtElement, except identity/generator constructors and predicates"
        if "fp12_in_group" in name or "in_group" in name:
            return "ReducedFp12"
        return "reduced Montgomery value of the declared field"
    if campaign == "groups":
        group = "G2" if any(token in name for token in ("p2", "g2")) else "G1"
        if "unchecked_mult" in name:
            return f"OnCurve{group} and ScalarBits"
        if "mult" in name:
            return f"Subgroup{group} and ScalarBits"
        if any(token in name for token in ("on_curve", "in_g1", "in_g2")):
            return f"FieldValid{group}Candidate"
        if re.search(r"(?:^|::|_)add(?:_affine)?$", name):
            return f"IncompleteAddPair<OnCurve{group}>"
        return f"OnCurve{group} in the declared representation"
    if campaign == "msm":
        return "bounded subgroup points, public ScalarBits, and typed disjoint scratch"
    if campaign == "pairing_state":
        if "miller" in name or "precompute_lines" in name:
            return "subgroup affine G1/G2 operands in the declared order"
        return "valid PairingState transition with subgroup operands"
    if campaign == "protocol":
        return "bounded byte slices or typed SecretScalar and subgroup operands"
    if campaign == "fft_middle":
        return "reduced Fr vectors, subgroup points, prepared domain, or live handle as declared"
    if campaign == "uniq":
        return "owned bounded messages and live uniqueness-tree storage"
    if campaign == "avx":
        return "typed scalar-equivalent value and supported AVX execution profile"
    return DOMAIN_BY_CAMPAIGN[campaign]


def factory_for(api: ApiSymbol, campaign: str, classification: str) -> str:
    name = api.symbol.lower()
    if classification == "indirect":
        return "corresponding_direct_surface"
    if campaign == "decode":
        return "bounded_decode_frame"
    if campaign == "bytes":
        return "bounded_bytes_words_or_typed_export"
    if campaign == "fields":
        if "final_exp" in name:
            return "nonzero_fp12_or_miller_factory"
        if "inverse" in name:
            return "nonzero_field_refinement"
        if "cyclotomic" in name or "gt_" in name or "::gt::" in name:
            return "gt_factory"
        return "reduced_field_factory"
    if campaign == "groups":
        if "unchecked_mult" in name:
            return "on_curve_point_factory"
        if "mult" in name:
            return "subgroup_point_factory"
        if re.search(r"(?:^|::|_)add(?:_affine)?$", name):
            return "incomplete_add_pair_factory"
        return "group_domain_factory"
    return FACTORY_BY_CAMPAIGN[campaign]


def row_for(api: ApiSymbol) -> dict[str, str]:
    campaign = campaign_for(api)
    classification = classification_for(api, campaign)
    if classification == "direct" and campaign in {"decode", "bytes"}:
        invalid = "malformed value or documented length rejection; readable storage remains valid"
    elif classification == "constructed" and campaign in {"msm", "fft_middle", "pairing_state", "uniq", "avx"}:
        invalid = "only documented size, scratch, overlap, or state rejection shapes"
    else:
        invalid = "none; factories skip unmet mathematical refinements"
    if api.surface == "rust_raw":
        notes = "raw Rust FFI inherits the corresponding C contract"
    elif api.surface == "cpp" and classification == "indirect":
        notes = "native method shares the implementation exercised by its direct owner surface"
    else:
        notes = "symbol disposition generated from the reviewed campaign and domain rules"
    return {
        "surface": api.surface,
        "symbol": api.symbol,
        "owner": api.owner,
        "campaign": campaign,
        "classification": classification,
        "required_domain": required_domain_for(api, campaign, classification),
        "factory": factory_for(api, campaign, classification),
        "allowed_invalid_cases": invalid,
        "oracles": ORACLE_BY_CAMPAIGN[campaign],
        "notes": notes,
    }


def load_config(root: Path) -> dict[str, Any]:
    return json.loads((root / "fuzz/api_inventory.json").read_text(encoding="utf-8"))


def collect_symbols(root: Path, require_generated: bool) -> list[ApiSymbol]:
    config = load_config(root)
    c_symbols = extract_c_symbols(root, config["c_headers"])
    symbols = list(c_symbols)
    symbols.extend(
        ApiSymbol("rust_raw", f"raw::{symbol.symbol}", symbol.owner)
        for symbol in c_symbols
    )

    cpp_ast = root / config["generated"]["cpp_ast"]
    rustdoc = root / config["generated"]["rustdoc"]
    if require_generated and not cpp_ast.is_file():
        raise ValueError(f"missing required Clang AST inventory: {cpp_ast}")
    if require_generated and not rustdoc.is_file():
        raise ValueError(f"missing required rustdoc inventory: {rustdoc}")
    cpp_snapshot = root / config["snapshots"]["cpp"]
    rust_snapshot = root / config["snapshots"]["rust_safe"]
    if cpp_ast.is_file():
        symbols.extend(
            extract_cpp_symbols(
                root,
                cpp_ast,
                config["cpp_owners"],
                config["cpp_public_types"],
            )
        )
    elif cpp_snapshot.is_file():
        symbols.extend(read_snapshot(cpp_snapshot, "cpp"))
    if rustdoc.is_file():
        symbols.extend(
            extract_rust_symbols(
                root,
                rustdoc,
                config["rust"]["manifest"],
                config["rust"]["excluded_modules"],
            )
        )
    elif rust_snapshot.is_file():
        symbols.extend(read_snapshot(rust_snapshot, "rust_safe"))
    return sorted(set(symbols))


def read_snapshot(path: Path, expected_surface: str) -> list[ApiSymbol]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema") != 1 or data.get("surface") != expected_surface:
        raise ValueError(f"invalid {expected_surface} API snapshot: {path}")
    result: list[ApiSymbol] = []
    for row in data.get("symbols", []):
        result.append(ApiSymbol(expected_surface, row["symbol"], row["owner"]))
    return result


def write_snapshot(path: Path, surface: str, symbols: Iterable[ApiSymbol]) -> None:
    selected = sorted(symbol for symbol in symbols if symbol.surface == surface)
    data = {
        "schema": 1,
        "surface": surface,
        "symbols": [
            {"symbol": symbol.symbol, "owner": symbol.owner}
            for symbol in selected
        ],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def write_ledger(path: Path, symbols: Iterable[ApiSymbol]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LEDGER_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for api in symbols:
            writer.writerow(row_for(api))


def read_ledger(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != LEDGER_COLUMNS:
            raise ValueError(
                f"ledger columns differ: expected {LEDGER_COLUMNS}, found {reader.fieldnames}"
            )
        return list(reader)


def validate_ledger(rows: list[dict[str, str]], expected: list[ApiSymbol]) -> list[str]:
    errors: list[str] = []
    keys: set[tuple[str, str]] = set()
    ledger_symbols: set[ApiSymbol] = set()
    for line_number, row in enumerate(rows, start=2):
        for column in LEDGER_COLUMNS:
            if not row[column].strip():
                errors.append(f"line {line_number}: empty {column}")
        key = (row["surface"], row["symbol"])
        if key in keys:
            errors.append(f"line {line_number}: duplicate symbol {key}")
        keys.add(key)
        if row["surface"] not in ALLOWED_SURFACES:
            errors.append(f"line {line_number}: unknown surface {row['surface']}")
        if row["campaign"] not in ALLOWED_CAMPAIGNS:
            errors.append(f"line {line_number}: unknown campaign {row['campaign']}")
        if row["classification"] not in ALLOWED_CLASSIFICATIONS:
            errors.append(
                f"line {line_number}: unknown classification {row['classification']}"
            )
        if row["classification"] == "excluded" and row["campaign"] != "excluded":
            errors.append(f"line {line_number}: excluded row must use excluded campaign")
        ledger_symbols.add(ApiSymbol(row["surface"], row["symbol"], row["owner"]))

    expected_set = set(expected)
    for missing in sorted(expected_set - ledger_symbols):
        errors.append(f"missing API row: {missing.surface}:{missing.symbol} ({missing.owner})")
    for extra in sorted(ledger_symbols - expected_set):
        errors.append(f"stale API row: {extra.surface}:{extra.symbol} ({extra.owner})")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument(
        "--ledger", type=Path, default=Path("docs/audit/FUZZ_API_LEDGER.csv")
    )
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--write-snapshots", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--require-generated", action="store_true")
    args = parser.parse_args()

    root = args.root.resolve()
    ledger = args.ledger if args.ledger.is_absolute() else root / args.ledger
    try:
        symbols = collect_symbols(root, args.require_generated)
        if args.self_test:
            sample = [row_for(symbols[0]), row_for(symbols[1])]
            missing_errors = validate_ledger(sample[:-1], symbols[:2])
            if not any("missing API row" in error for error in missing_errors):
                raise ValueError("self-test did not reject a missing API row")
            duplicate = [dict(sample[0]), dict(sample[0])]
            duplicate_errors = validate_ledger(duplicate, symbols[:1])
            if not any("duplicate symbol" in error for error in duplicate_errors):
                raise ValueError("self-test did not reject a duplicate API row")
            unknown = [dict(sample[0])]
            unknown[0]["campaign"] = "unknown"
            unknown_errors = validate_ledger(unknown, symbols[:1])
            if not any("unknown campaign" in error for error in unknown_errors):
                raise ValueError("self-test did not reject an unknown campaign")
        if args.write_snapshots:
            config = load_config(root)
            for surface in ("cpp", "rust_safe"):
                if any(symbol.surface == surface for symbol in symbols):
                    write_snapshot(root / config["snapshots"][surface], surface, symbols)
            print("updated generated API snapshots")
        if args.write:
            write_ledger(ledger, symbols)
            print(f"wrote {len(symbols)} API rows to {ledger}")
            return 0
        rows = read_ledger(ledger)
        errors = validate_ledger(rows, symbols)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"fuzz API inventory failed: {error}", file=sys.stderr)
        return 1

    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        print(f"fuzz API inventory failed with {len(errors)} error(s)", file=sys.stderr)
        return 1
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["surface"]] = counts.get(row["surface"], 0) + 1
    summary = ", ".join(f"{surface}={count}" for surface, count in sorted(counts.items()))
    print(f"fuzz API ledger passed: {len(rows)} rows ({summary})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
