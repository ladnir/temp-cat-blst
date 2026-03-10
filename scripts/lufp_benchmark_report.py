#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import pathlib
import subprocess
import tempfile
from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Tuple


CANONICAL_ADD_BASELINE = {
    "BM_LUFpAddAccounting<LUFpAddStage::RawAdd>/256": 12216.0,
    "BM_LUFpAddAccounting<LUFpAddStage::RawAddBias>/256": 19346.0,
    "BM_LUFpAddAccounting<LUFpAddStage::RawAddBiasCarry>/256": 49990.0,
    "BM_LUFpAddAccounting<LUFpAddStage::Full>/256": 81305.0,
    "BM_LUFpAddAccounting<LUFpAddStage::FullSign4>/256": 181537.0,
}

REQUIRED_REPETITIONS = 7
CANONICAL_ADD_TOLERANCE = 0.03

FAMILIES = {
    "UFp add accounting": [
        "BM_UFpAddAccounting<UFpAddStage::RawAdd>/256",
        "BM_UFpAddAccounting<UFpAddStage::RawAddBias>/256",
        "BM_UFpAddAccounting<UFpAddStage::RawAddBiasCarry>/256",
        "BM_UFpAddAccounting<UFpAddStage::Full>/256",
    ],
    "LUFp v2 add accounting": [
        "BM_LUFpAddAccounting<LUFpAddStage::RawAdd>/256",
        "BM_LUFpAddAccounting<LUFpAddStage::RawAddBias>/256",
        "BM_LUFpAddAccounting<LUFpAddStage::RawAddBiasCarry>/256",
        "BM_LUFpAddAccounting<LUFpAddStage::Full>/256",
        "BM_LUFpAddAccounting<LUFpAddStage::FullSign4>/256",
    ],
    "LUFp v2 renorm accounting": [
        "BM_LUFpV2ExactAddRenormAccounting<LUFpV2ExactAddRenormStage::Bias>/256",
        "BM_LUFpV2ExactAddRenormAccounting<LUFpV2ExactAddRenormStage::BiasCarry>/256",
        "BM_LUFpV2ExactAddRenormAccounting<LUFpV2ExactAddRenormStage::Full>/256",
        "BM_LUFpV2ExactAddRenormAccounting<LUFpV2ExactAddRenormStage::FullSign4>/256",
    ],
    "UFp Fp6 mul accounting": [
        "BM_Fp6MulUfpAccounting<UFpFp6MulStage::PreAdd>/256",
        "BM_Fp6MulUfpAccounting<UFpFp6MulStage::PreAddMul>/256",
        "BM_Fp6MulUfpAccounting<UFpFp6MulStage::Full>/256",
    ],
    "LUFp v2 Fp6 mul accounting": [
        "BM_Fp6MulLufpAccounting<LUFpFp6MulStage::PreAdd>/256",
        "BM_Fp6MulLufpAccounting<LUFpFp6MulStage::PreAddMul>/256",
        "BM_Fp6MulLufpAccounting<LUFpFp6MulStage::PreAddMulRawCombine>/256",
        "BM_Fp6MulLufpAccounting<LUFpFp6MulStage::Full>/256",
    ],
    "LUFp v2 Fp6 renorm accounting": [
        "BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::CopyBias>/256",
        "BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::CopyBiasCarry>/256",
        "BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::CopyBiasCarrySign4>/256",
        "BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::CopyBiasCarrySign42>/256",
        "BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::CopyBiasCarrySign421>/256",
        "BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::CopyBiasCarrySign4211>/256",
        "BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::Full>/256",
    ],
}

PROFILES = {
    "canonical-add-only": {
        "required": [
            "UFp add accounting",
            "LUFp v2 add accounting",
        ],
        "filter": r"^(BM_UFpAddAccounting.*|BM_LUFpAddAccounting.*)$",
    },
    "canonical-core": {
        "required": [
            "UFp add accounting",
            "LUFp v2 add accounting",
            "LUFp v2 renorm accounting",
        ],
        "filter": r"^(BM_UFpAddAccounting.*|BM_LUFpAddAccounting.*|BM_LUFpV2ExactAddRenormAccounting.*)$",
    },
    "canonical-full": {
        "required": [
            "UFp add accounting",
            "LUFp v2 add accounting",
            "LUFp v2 renorm accounting",
            "UFp Fp6 mul accounting",
            "LUFp v2 Fp6 mul accounting",
        ],
        "filter": (
            r"^(BM_UFpAddAccounting.*|"
            r"BM_LUFpAddAccounting.*|"
            r"BM_LUFpV2ExactAddRenormAccounting.*|"
            r"BM_Fp6MulUfpAccounting.*|"
            r"BM_Fp6MulLufpAccounting.*|"
            r"BM_LUFpV2Fp6MulCommonRenormAccounting.*)$"
        ),
    },
    "exploratory": {
        "required": [],
        "filter": None,
    },
}


@dataclass
class Aggregate:
    mean: float
    median: float
    stddev: float
    cv: float
    time_unit: str


def run(cmd: List[str], cwd: pathlib.Path) -> str:
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"command failed: {' '.join(cmd)}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    return proc.stdout.strip()


def git_metadata(repo: pathlib.Path) -> Tuple[str, str]:
    commit = run(["git", "rev-parse", "HEAD"], repo)
    dirty = "dirty" if run(["git", "status", "--porcelain"], repo) else "clean"
    return commit, dirty


def cpu_model() -> str:
    cpuinfo = pathlib.Path("/proc/cpuinfo")
    if cpuinfo.exists():
        for line in cpuinfo.read_text().splitlines():
            if line.lower().startswith("model name"):
                return line.split(":", 1)[1].strip()
    return "unknown"


def load_average() -> Tuple[float, float, float]:
    try:
        values = os.getloadavg()
        return float(values[0]), float(values[1]), float(values[2])
    except OSError:
        return 0.0, 0.0, 0.0


def parse_benchmark_json(path: pathlib.Path) -> Dict[str, Aggregate]:
    payload = json.loads(path.read_text())
    grouped: Dict[str, Dict[str, Tuple[float, str]]] = {}
    for entry in payload.get("benchmarks", []):
        agg = entry.get("aggregate_name")
        if agg not in {"mean", "median", "stddev"}:
            continue
        name = entry["name"]
        suffix = f"_{agg}"
        base = name[:-len(suffix)] if name.endswith(suffix) else name
        grouped.setdefault(base, {})[agg] = (float(entry["real_time"]), entry["time_unit"])

    out: Dict[str, Aggregate] = {}
    for base, fields in grouped.items():
        if {"mean", "median", "stddev"} - set(fields):
            continue
        mean, unit = fields["mean"]
        median, _ = fields["median"]
        stddev, _ = fields["stddev"]
        cv = 0.0 if mean == 0.0 else (stddev / mean) * 100.0
        out[base] = Aggregate(mean=mean, median=median, stddev=stddev, cv=cv, time_unit=unit)
    return out


def fmt_ns(value: float) -> str:
    return f"{value:.2f} ns"


def fmt_pct(value: float) -> str:
    return f"{value:.2f}%"


def add_table(lines: List[str], headers: List[str], rows: List[List[str]]) -> None:
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("|" + "|".join(["---"] * len(headers)) + "|")
    for row in rows:
        lines.append("| " + " | ".join(row) + " |")
    lines.append("")


def family_present(results: Dict[str, Aggregate], family: str) -> bool:
    return all(name in results for name in FAMILIES[family])


def stage_rows(results: Dict[str, Aggregate], names: List[str]) -> List[Tuple[str, float, float, float]]:
    rows: List[Tuple[str, float, float, float]] = []
    prev = 0.0
    full = results[names[-1]].mean
    for name in names:
        cumulative = results[name].mean
        delta = cumulative - prev
        pct = 0.0 if full == 0.0 else (delta / full) * 100.0
        rows.append((name, cumulative, delta, pct))
        prev = cumulative
    return rows


def alt_terminal_rows(
    results: Dict[str, Aggregate],
    names: List[str],
    alternative_terminal: str,
) -> List[Tuple[str, float, float, float]]:
    rows: List[Tuple[str, float, float, float]] = []
    prev = 0.0
    full = results[alternative_terminal].mean
    for name in names:
        cumulative = results[name].mean
        delta = cumulative - prev
        pct = 0.0 if full == 0.0 else (delta / full) * 100.0
        rows.append((name, cumulative, delta, pct))
        prev = cumulative
    predecessor = results[names[-1]].mean
    delta = full - predecessor
    pct = 0.0 if full == 0.0 else (delta / full) * 100.0
    rows.append((alternative_terminal, full, delta, pct))
    return rows


def canonical_add_flags(results: Dict[str, Aggregate]) -> List[Tuple[str, str, str, str]]:
    flags: List[Tuple[str, str, str, str]] = []
    if not family_present(results, "LUFp v2 add accounting"):
        return flags

    for name, baseline in CANONICAL_ADD_BASELINE.items():
        measured = results[name].mean
        if abs(measured - baseline) / baseline > CANONICAL_ADD_TOLERANCE:
            flags.append(
                (
                    "fatal",
                    "LUFp v2 add accounting",
                    f"{name} measured {fmt_ns(measured)} vs baseline {fmt_ns(baseline)}",
                    "revert the regression or explain and accept a new baseline explicitly",
                )
            )

    raw = results["BM_LUFpAddAccounting<LUFpAddStage::RawAdd>/256"].mean
    bias = results["BM_LUFpAddAccounting<LUFpAddStage::RawAddBias>/256"].mean
    carry = results["BM_LUFpAddAccounting<LUFpAddStage::RawAddBiasCarry>/256"].mean
    full = results["BM_LUFpAddAccounting<LUFpAddStage::Full>/256"].mean
    full_sign4 = results["BM_LUFpAddAccounting<LUFpAddStage::FullSign4>/256"].mean
    add_ops = 256.0 * 16.0

    canaries = {
        "bias delta": (bias - raw, 1.74),
        "carry delta": (carry - bias, 7.48),
        "1-step repair delta": (full - carry, 7.65),
        "4-step repair delta": (full_sign4 - carry, 32.12),
    }
    for label, (measured, baseline) in canaries.items():
        measured_per_add = measured / add_ops
        if abs(measured_per_add - baseline) / baseline > CANONICAL_ADD_TOLERANCE:
            flags.append(
                (
                    "fatal",
                    "LUFp v2 add accounting",
                    f"{label} measured {measured_per_add:.2f} ns/add vs baseline {baseline:.2f} ns/add",
                    "revert the regression or accept a new baseline explicitly",
                )
            )

    return flags


def collect_red_flags(
    results: Dict[str, Aggregate],
    profile: str,
    idle_threshold: float,
    load1: float,
    repetitions: int,
    min_time: str,
) -> List[Tuple[str, str, str, str]]:
    flags: List[Tuple[str, str, str, str]] = []

    for family in PROFILES[profile]["required"]:
        missing = [name for name in FAMILIES[family] if name not in results]
        if missing:
            flags.append(
                (
                    "fatal",
                    family,
                    f"missing benchmarks: {', '.join(missing)}",
                    "fix benchmark registration or the selected filter before trusting the run",
                )
            )

    for name, agg in results.items():
        if agg.cv > 5.0:
            flags.append(("fatal", name, f"coefficient of variation is {fmt_pct(agg.cv)}", "rerun on a more stable idle machine"))
        elif agg.cv > 2.0:
            flags.append(("warning", name, f"coefficient of variation is {fmt_pct(agg.cv)}", "rerun to improve stability"))

    if repetitions < REQUIRED_REPETITIONS:
        flags.append(("warning", "run", f"repetition count {repetitions} is below the preferred {REQUIRED_REPETITIONS}", "use the standard repetition count for canonical runs"))

    if min_time in {"0", "0s", "0.0", "0.0s"}:
        flags.append(("warning", "run", "minimum benchmark time is effectively zero", "use a non-trivial benchmark min time"))

    if load1 > idle_threshold * 2.0:
        flags.append(("fatal", "environment", f"load average before run was {load1:.2f} (threshold {idle_threshold:.2f})", "rerun on an idle machine"))
    elif load1 > idle_threshold:
        flags.append(("warning", "environment", f"load average before run was {load1:.2f} (threshold {idle_threshold:.2f})", "prefer an idler machine or rerun"))

    if profile != "exploratory":
        flags.extend(canonical_add_flags(results))

    return flags


def speed_of_light_lines(results: Dict[str, Aggregate]) -> List[str]:
    lines: List[str] = []
    families = [
        (
            "UFp add accounting",
            "BM_UFpAddAccounting<UFpAddStage::RawAdd>/256",
            "BM_UFpAddAccounting<UFpAddStage::Full>/256",
            "Lower bound is the measured raw-add stage. Assumption: later stages are non-negative overhead.",
        ),
        (
            "LUFp v2 add accounting",
            "BM_LUFpAddAccounting<LUFpAddStage::RawAdd>/256",
            "BM_LUFpAddAccounting<LUFpAddStage::Full>/256",
            "Lower bound is the measured raw lazy add stage. Assumption: bias, carry, and repair are non-negative overhead.",
        ),
        (
            "LUFp v2 renorm accounting",
            "BM_LUFpV2ExactAddRenormAccounting<LUFpV2ExactAddRenormStage::Bias>/256",
            "BM_LUFpV2ExactAddRenormAccounting<LUFpV2ExactAddRenormStage::Full>/256",
            "Lower bound is the measured bias stage. Assumption: carry and repair are non-negative overhead.",
        ),
        (
            "UFp Fp6 mul accounting",
            "BM_Fp6MulUfpAccounting<UFpFp6MulStage::PreAddMul>/256",
            "BM_Fp6MulUfpAccounting<UFpFp6MulStage::Full>/256",
            "Lower bound is the measured pre-add plus multiply stage. Assumption: the remaining combine/output work is non-negative overhead.",
        ),
        (
            "LUFp v2 Fp6 mul accounting",
            "BM_Fp6MulLufpAccounting<LUFpFp6MulStage::PreAddMulRawCombine>/256",
            "BM_Fp6MulLufpAccounting<LUFpFp6MulStage::Full>/256",
            "Lower bound is the measured pre-add plus multiply plus raw-combine stage. Assumption: the final boundary is non-negative overhead.",
        ),
        (
            "LUFp v2 Fp6 renorm accounting",
            "BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::CopyBias>/256",
            "BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::Full>/256",
            "Lower bound is the measured copy-plus-bias stage. Assumption: carry, repair, and output are non-negative overhead.",
        ),
    ]

    for label, lower_name, measured_name, assumptions in families:
        if lower_name not in results or measured_name not in results:
            continue
        lower = results[lower_name].mean
        measured = results[measured_name].mean
        ratio = 0.0 if lower == 0.0 else measured / lower
        gap = measured - lower
        lines.append(f"### {label}")
        lines.append(f"- Lower bound: {fmt_ns(lower)}")
        lines.append(f"- Measured: {fmt_ns(measured)}")
        lines.append(f"- Gap: {fmt_ns(gap)} ({ratio:.2f}x of lower bound)")
        lines.append(f"- Assumptions: {assumptions}")
        lines.append("- Facts: lower and measured values come directly from this run.")
        lines.append("- Inferences: the gap is attributable to stages excluded from the chosen lower-bound stage.")
        lines.append("- Speculation: none beyond the stated stage decomposition.")
        lines.append("")
    return lines


def per_operation_rows(results: Dict[str, Aggregate]) -> List[List[str]]:
    rows: List[List[str]] = []
    add_ops = 256 * 16

    if family_present(results, "UFp add accounting"):
        ufp_full = results["BM_UFpAddAccounting<UFpAddStage::Full>/256"].mean
        rows.append(["UFp add", "full strict add", fmt_ns(ufp_full / add_ops)])

    if family_present(results, "LUFp v2 add accounting"):
        raw = results["BM_LUFpAddAccounting<LUFpAddStage::RawAdd>/256"].mean
        bias = results["BM_LUFpAddAccounting<LUFpAddStage::RawAddBias>/256"].mean
        carry = results["BM_LUFpAddAccounting<LUFpAddStage::RawAddBiasCarry>/256"].mean
        full = results["BM_LUFpAddAccounting<LUFpAddStage::Full>/256"].mean
        full_sign4 = results["BM_LUFpAddAccounting<LUFpAddStage::FullSign4>/256"].mean
        rows.append(["LUFp v2 add", "raw add", fmt_ns(raw / add_ops)])
        rows.append(["LUFp v2 add", "bias delta", fmt_ns((bias - raw) / add_ops)])
        rows.append(["LUFp v2 add", "carry delta", fmt_ns((carry - bias) / add_ops)])
        rows.append(["LUFp v2 add", "1-step repair delta", fmt_ns((full - carry) / add_ops)])
        rows.append(["LUFp v2 add", "4-step repair delta", fmt_ns((full_sign4 - carry) / add_ops)])

    if family_present(results, "LUFp v2 renorm accounting"):
        bias = results["BM_LUFpV2ExactAddRenormAccounting<LUFpV2ExactAddRenormStage::Bias>/256"].mean
        carry = results["BM_LUFpV2ExactAddRenormAccounting<LUFpV2ExactAddRenormStage::BiasCarry>/256"].mean
        full = results["BM_LUFpV2ExactAddRenormAccounting<LUFpV2ExactAddRenormStage::Full>/256"].mean
        full_sign4 = results["BM_LUFpV2ExactAddRenormAccounting<LUFpV2ExactAddRenormStage::FullSign4>/256"].mean
        rows.append(["LUFp v2 renorm", "bias", fmt_ns(bias / add_ops)])
        rows.append(["LUFp v2 renorm", "carry delta", fmt_ns((carry - bias) / add_ops)])
        rows.append(["LUFp v2 renorm", "1-step repair delta", fmt_ns((full - carry) / add_ops)])
        rows.append(["LUFp v2 renorm", "4-step repair delta", fmt_ns((full_sign4 - carry) / add_ops)])

    if family_present(results, "UFp Fp6 mul accounting"):
        pre_add = results["BM_Fp6MulUfpAccounting<UFpFp6MulStage::PreAdd>/256"].mean
        pre_add_mul = results["BM_Fp6MulUfpAccounting<UFpFp6MulStage::PreAddMul>/256"].mean
        full = results["BM_Fp6MulUfpAccounting<UFpFp6MulStage::Full>/256"].mean
        rows.append(["UFp Fp6 mul", "narrow pre-add (6 Fp2 adds)", fmt_ns(pre_add / (256 * 6))])
        rows.append(["UFp Fp6 mul", "Fp2 mul (6 muls)", fmt_ns((pre_add_mul - pre_add) / (256 * 6))])
        rows.append(["UFp Fp6 mul", "strict combine/output boundary", fmt_ns((full - pre_add_mul) / 256)])
        rows.append(["UFp Fp6 mul", "combine add/sub-like op (13 Fp2 ops)", fmt_ns((full - pre_add_mul) / (256 * 13))])

    if family_present(results, "LUFp v2 Fp6 mul accounting"):
        pre_add = results["BM_Fp6MulLufpAccounting<LUFpFp6MulStage::PreAdd>/256"].mean
        pre_add_mul = results["BM_Fp6MulLufpAccounting<LUFpFp6MulStage::PreAddMul>/256"].mean
        raw_combine = results["BM_Fp6MulLufpAccounting<LUFpFp6MulStage::PreAddMulRawCombine>/256"].mean
        full = results["BM_Fp6MulLufpAccounting<LUFpFp6MulStage::Full>/256"].mean
        rows.append(["LUFp v2 Fp6 mul", "narrow pre-add (6 Fp2 adds)", fmt_ns(pre_add / (256 * 6))])
        rows.append(["LUFp v2 Fp6 mul", "Fp2 mul (6 muls)", fmt_ns((pre_add_mul - pre_add) / (256 * 6))])
        rows.append(["LUFp v2 Fp6 mul", "raw lazy combine", fmt_ns((raw_combine - pre_add_mul) / 256)])
        rows.append(["LUFp v2 Fp6 mul", "raw combine op (13 LUFp2 ops)", fmt_ns((raw_combine - pre_add_mul) / (256 * 13))])
        rows.append(["LUFp v2 Fp6 mul", "final renorm/output boundary", fmt_ns((full - raw_combine) / 256)])

    if family_present(results, "LUFp v2 Fp6 renorm accounting"):
        bias = results["BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::CopyBias>/256"].mean
        carry = results["BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::CopyBiasCarry>/256"].mean
        sign4 = results["BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::CopyBiasCarrySign4>/256"].mean
        sign42 = results["BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::CopyBiasCarrySign42>/256"].mean
        sign421 = results["BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::CopyBiasCarrySign421>/256"].mean
        sign4211 = results["BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::CopyBiasCarrySign4211>/256"].mean
        full = results["BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::Full>/256"].mean
        rows.append(["LUFp v2 Fp6 renorm", "copy+bias (6 components)", fmt_ns(bias / (256 * 6))])
        rows.append(["LUFp v2 Fp6 renorm", "carry delta (6 components)", fmt_ns((carry - bias) / (256 * 6))])
        rows.append(["LUFp v2 Fp6 renorm", "sign4 delta (6 components)", fmt_ns((sign4 - carry) / (256 * 6))])
        rows.append(["LUFp v2 Fp6 renorm", "sign2 delta (6 components)", fmt_ns((sign42 - sign4) / (256 * 6))])
        rows.append(["LUFp v2 Fp6 renorm", "sign1 delta A (6 components)", fmt_ns((sign421 - sign42) / (256 * 6))])
        rows.append(["LUFp v2 Fp6 renorm", "sign1 delta B (6 components)", fmt_ns((sign4211 - sign421) / (256 * 6))])
        rows.append(["LUFp v2 Fp6 renorm", "store/output delta (6 components)", fmt_ns((full - sign4211) / (256 * 6))])

    return rows


def fp6_operation_accounting_lines(results: Dict[str, Aggregate]) -> List[str]:
    lines: List[str] = []

    if family_present(results, "UFp Fp6 mul accounting"):
        pre_add = results["BM_Fp6MulUfpAccounting<UFpFp6MulStage::PreAdd>/256"].mean / 256.0
        pre_add_mul = results["BM_Fp6MulUfpAccounting<UFpFp6MulStage::PreAddMul>/256"].mean / 256.0
        full = results["BM_Fp6MulUfpAccounting<UFpFp6MulStage::Full>/256"].mean / 256.0
        mul_total = pre_add_mul - pre_add
        combine_total = full - pre_add_mul
        lines.append("### UFp Fp6 operation model")
        add_table(
            lines,
            ["Bucket", "Count per Fp6", "Total per Fp6", "Average per counted op", "Notes"],
            [
                ["Narrow pre-add", "6 Fp2 adds", fmt_ns(pre_add), fmt_ns(pre_add / 6.0), "Measured by `PreAdd`."],
                ["Fp2 mul", "6 Fp2 muls", fmt_ns(mul_total), fmt_ns(mul_total / 6.0), "Measured by `PreAddMul - PreAdd`."],
                ["Strict combine/output", "13 Fp2 add/sub-like ops", fmt_ns(combine_total), fmt_ns(combine_total / 13.0), "Derived from the current Karatsuba recombination shape in the benchmark code."],
                ["Model sum", "1 Fp6 mul", fmt_ns(pre_add + mul_total + combine_total), fmt_ns(pre_add + mul_total + combine_total), "This decomposition is stage-derived, so it closes exactly to the measured full stage."],
                ["Measured full", "1 Fp6 mul", fmt_ns(full), fmt_ns(full), "Measured by `Full`."],
            ],
        )

    if family_present(results, "LUFp v2 Fp6 mul accounting") and family_present(results, "LUFp v2 Fp6 renorm accounting"):
        pre_add = results["BM_Fp6MulLufpAccounting<LUFpFp6MulStage::PreAdd>/256"].mean / 256.0
        pre_add_mul = results["BM_Fp6MulLufpAccounting<LUFpFp6MulStage::PreAddMul>/256"].mean / 256.0
        raw_combine = results["BM_Fp6MulLufpAccounting<LUFpFp6MulStage::PreAddMulRawCombine>/256"].mean / 256.0
        full = results["BM_Fp6MulLufpAccounting<LUFpFp6MulStage::Full>/256"].mean / 256.0
        renorm_bias = results["BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::CopyBias>/256"].mean / 256.0
        renorm_carry = results["BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::CopyBiasCarry>/256"].mean / 256.0
        renorm_sign4 = results["BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::CopyBiasCarrySign4>/256"].mean / 256.0
        renorm_sign42 = results["BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::CopyBiasCarrySign42>/256"].mean / 256.0
        renorm_sign421 = results["BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::CopyBiasCarrySign421>/256"].mean / 256.0
        renorm_sign4211 = results["BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::CopyBiasCarrySign4211>/256"].mean / 256.0
        renorm_full = results["BM_LUFpV2Fp6MulCommonRenormAccounting<LUFpV2Fp6MulCommonRenormStage::Full>/256"].mean / 256.0

        mul_total = pre_add_mul - pre_add
        raw_total = raw_combine - pre_add_mul
        bias_total = renorm_bias
        carry_total = renorm_carry - renorm_bias
        sign4_total = renorm_sign4 - renorm_carry
        sign2_total = renorm_sign42 - renorm_sign4
        sign1a_total = renorm_sign421 - renorm_sign42
        sign1b_total = renorm_sign4211 - renorm_sign421
        store_total = renorm_full - renorm_sign4211
        model_sum = pre_add + mul_total + raw_total + renorm_full
        delta = model_sum - full
        pct = 0.0 if full == 0.0 else (delta / full) * 100.0

        lines.append("### LUFp v2 Fp6 operation model")
        add_table(
            lines,
            ["Bucket", "Count per Fp6", "Total per Fp6", "Average per counted op", "Notes"],
            [
                ["Narrow pre-add", "6 Fp2 adds", fmt_ns(pre_add), fmt_ns(pre_add / 6.0), "Measured by `PreAdd`."],
                ["Fp2 mul", "6 Fp2 muls", fmt_ns(mul_total), fmt_ns(mul_total / 6.0), "Measured by `PreAddMul - PreAdd`."],
                ["Raw lazy combine", "13 LUFp2 add/sub-like ops", fmt_ns(raw_total), fmt_ns(raw_total / 13.0), "Measured by `PreAddMulRawCombine - PreAddMul`."],
                ["Renorm copy+bias", "6 components", fmt_ns(bias_total), fmt_ns(bias_total / 6.0), "Measured independently on raw-combine outputs."],
                ["Renorm carry", "6 components", fmt_ns(carry_total), fmt_ns(carry_total / 6.0), "Measured independently on raw-combine outputs."],
                ["Renorm sign4", "6 components", fmt_ns(sign4_total), fmt_ns(sign4_total / 6.0), "Measured independently on raw-combine outputs."],
                ["Renorm sign2", "6 components", fmt_ns(sign2_total), fmt_ns(sign2_total / 6.0), "Measured independently on raw-combine outputs."],
                ["Renorm sign1 A", "6 components", fmt_ns(sign1a_total), fmt_ns(sign1a_total / 6.0), "Measured independently on raw-combine outputs."],
                ["Renorm sign1 B", "6 components", fmt_ns(sign1b_total), fmt_ns(sign1b_total / 6.0), "Measured independently on raw-combine outputs."],
                ["Renorm store/output", "6 components", fmt_ns(store_total), fmt_ns(store_total / 6.0), "Measured independently on raw-combine outputs."],
                ["Model sum", "1 Fp6 mul", fmt_ns(model_sum), fmt_ns(model_sum), "Uses independently measured renorm buckets plus the main Fp6 stages."],
                ["Measured full", "1 Fp6 mul", fmt_ns(full), fmt_ns(full), "Measured by `BM_Fp6MulLufpAccounting<...Full>`."],
                ["Model error", "relative to measured full", fmt_ns(delta), fmt_pct(pct), "Positive means the model overestimates the measured full path."],
            ],
        )

    return lines


def render_report(
    repo: pathlib.Path,
    results: Dict[str, Aggregate],
    output_path: pathlib.Path,
    command: str,
    profile: str,
    repetitions: int,
    min_time: str,
    idle_threshold: float,
    server_ip: str,
    affinity: str,
) -> None:
    timestamp = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")
    commit, dirty = git_metadata(repo)
    load1, load5, load15 = load_average()
    flags = collect_red_flags(results, profile, idle_threshold, load1, repetitions, min_time)

    lines: List[str] = []
    lines.append("# LUFp Benchmark Report")
    lines.append("")
    lines.append("## Red Flags")
    if not flags:
        lines.append("None")
    else:
        for severity, family, explanation, action in flags:
            lines.append(f"- {severity}: `{family}`: {explanation}. Suggested action: {action}.")
    lines.append("")

    lines.append("## Run Metadata")
    lines.append(f"- Timestamp: {timestamp}")
    lines.append(f"- Run profile: `{profile}`")
    lines.append(f"- Git commit: `{commit}`")
    lines.append(f"- Git status: `{dirty}`")
    lines.append(f"- Hostname: `{os.uname().nodename}`")
    if server_ip:
        lines.append(f"- Server IP: `{server_ip}`")
    lines.append(f"- CPU model: {cpu_model()}")
    lines.append(f"- Load average before run: {load1:.2f}, {load5:.2f}, {load15:.2f}")
    lines.append(f"- Benchmark binary path: `{command.split()[0]}`")
    lines.append(f"- Benchmark command line: `{command}`")
    lines.append(f"- Repetition count: `{repetitions}`")
    lines.append(f"- Minimum benchmark time: `{min_time}`")
    lines.append(f"- CPU affinity: `{affinity or 'unspecified'}`")
    lines.append("")

    lines.append("## Benchmark Coverage")
    required = PROFILES[profile]["required"]
    for family in required:
        missing = [name for name in FAMILIES[family] if name not in results]
        status = "present" if not missing else "missing"
        lines.append(f"- {family}: {status}")
        if missing:
            lines.append(f"  missing benchmarks: {', '.join(missing)}")
    if profile == "exploratory":
        lines.append("- exploratory: report only the families intentionally targeted by the selected filter")
    lines.append("")

    lines.append("## Raw Results")
    raw_rows = []
    for name in sorted(results.keys()):
        agg = results[name]
        raw_rows.append([name, fmt_ns(agg.mean), fmt_ns(agg.median), fmt_ns(agg.stddev), fmt_pct(agg.cv)])
    add_table(lines, ["Benchmark", "Mean", "Median", "Stddev", "CV"], raw_rows)

    lines.append("## Speed-of-Light Summary")
    lines.extend(speed_of_light_lines(results))

    lines.append("## Accounting Breakdown")
    if family_present(results, "UFp add accounting"):
        lines.append("### UFp add accounting")
        add_table(
            lines,
            ["Stage", "Cumulative", "Delta", "% Full"],
            [[name, fmt_ns(c), fmt_ns(d), fmt_pct(p)] for name, c, d, p in stage_rows(results, FAMILIES["UFp add accounting"])],
        )

    if family_present(results, "LUFp v2 add accounting"):
        exact = FAMILIES["LUFp v2 add accounting"][:-1]
        alt = FAMILIES["LUFp v2 add accounting"][-1]
        lines.append("### LUFp v2 add accounting (ExactAdd)")
        add_table(
            lines,
            ["Stage", "Cumulative", "Delta", "% Full"],
            [[name, fmt_ns(c), fmt_ns(d), fmt_pct(p)] for name, c, d, p in stage_rows(results, exact)],
        )
        lines.append("### LUFp v2 add accounting (FullSign4)")
        add_table(
            lines,
            ["Stage", "Cumulative", "Derived Delta", "% FullSign4"],
            [[name, fmt_ns(c), fmt_ns(d), fmt_pct(p)] for name, c, d, p in alt_terminal_rows(results, exact[:-1], alt)],
        )

    if family_present(results, "LUFp v2 renorm accounting"):
        exact = FAMILIES["LUFp v2 renorm accounting"][:-1]
        alt = FAMILIES["LUFp v2 renorm accounting"][-1]
        lines.append("### LUFp v2 renorm accounting (ExactAdd)")
        add_table(
            lines,
            ["Stage", "Cumulative", "Delta", "% Full"],
            [[name, fmt_ns(c), fmt_ns(d), fmt_pct(p)] for name, c, d, p in stage_rows(results, exact)],
        )
        lines.append("### LUFp v2 renorm accounting (FullSign4)")
        add_table(
            lines,
            ["Stage", "Cumulative", "Derived Delta", "% FullSign4"],
            [[name, fmt_ns(c), fmt_ns(d), fmt_pct(p)] for name, c, d, p in alt_terminal_rows(results, exact[:-1], alt)],
        )

    if family_present(results, "UFp Fp6 mul accounting"):
        lines.append("### UFp Fp6 mul accounting")
        add_table(
            lines,
            ["Stage", "Cumulative", "Delta", "% Full"],
            [[name, fmt_ns(c), fmt_ns(d), fmt_pct(p)] for name, c, d, p in stage_rows(results, FAMILIES["UFp Fp6 mul accounting"])],
        )

    if family_present(results, "LUFp v2 Fp6 mul accounting"):
        lines.append("### LUFp v2 Fp6 mul accounting")
        add_table(
            lines,
            ["Stage", "Cumulative", "Delta", "% Full"],
            [[name, fmt_ns(c), fmt_ns(d), fmt_pct(p)] for name, c, d, p in stage_rows(results, FAMILIES["LUFp v2 Fp6 mul accounting"])],
        )

    if family_present(results, "LUFp v2 Fp6 renorm accounting"):
        lines.append("### LUFp v2 Fp6 renorm accounting")
        add_table(
            lines,
            ["Stage", "Cumulative", "Delta", "% Full"],
            [[name, fmt_ns(c), fmt_ns(d), fmt_pct(p)] for name, c, d, p in stage_rows(results, FAMILIES["LUFp v2 Fp6 renorm accounting"])],
        )

    lines.append("## Per-Operation Summary")
    per_rows = per_operation_rows(results)
    if per_rows:
        add_table(lines, ["Family", "Operation", "Average"], per_rows)
    else:
        lines.append("None")
        lines.append("")

    lines.append("## Fp6 Operation Accounting")
    fp6_lines = fp6_operation_accounting_lines(results)
    if fp6_lines:
        lines.extend(fp6_lines)
    else:
        lines.append("None")
        lines.append("")

    lines.append("## Full Comparisons")
    full_rows: List[List[str]] = []
    if family_present(results, "UFp add accounting") and family_present(results, "LUFp v2 add accounting"):
        a = results["BM_UFpAddAccounting<UFpAddStage::Full>/256"].mean
        b = results["BM_LUFpAddAccounting<LUFpAddStage::Full>/256"].mean
        full_rows.append(["UFp add vs LUFp v2 add (Full)", fmt_ns(a), fmt_ns(b), fmt_ns(b - a), fmt_pct(((b - a) / a) * 100.0)])
    if family_present(results, "UFp Fp6 mul accounting") and family_present(results, "LUFp v2 Fp6 mul accounting"):
        a = results["BM_Fp6MulUfpAccounting<UFpFp6MulStage::Full>/256"].mean
        b = results["BM_Fp6MulLufpAccounting<LUFpFp6MulStage::Full>/256"].mean
        full_rows.append(["UFp Fp6 mul vs LUFp v2 Fp6 mul (Full)", fmt_ns(a), fmt_ns(b), fmt_ns(b - a), fmt_pct(((b - a) / a) * 100.0)])
    if full_rows:
        add_table(lines, ["Comparison", "Baseline", "Candidate", "Absolute Delta", "Percent Delta"], full_rows)
    else:
        lines.append("None")
        lines.append("")

    lines.append("## Interpretation")
    if any(severity == "fatal" for severity, _, _, _ in flags):
        lines.append("- Run not trustworthy for final conclusions")
    lines.append("- Facts: this report treats accounting benchmarks as the source of constituent timing.")
    lines.append("- Inferences: any lower-bound gap reported above is attributable only to stages excluded from the chosen lower-bound stage.")
    lines.append("- Speculation: none beyond the explicitly labeled speed-of-light sections.")
    lines.append("")

    lines.append("## Notes")
    lines.append("- Canonical add-family baselines are enforced for non-exploratory profiles.")
    lines.append("- Alternative `FullSign4` repair costs are reported relative to the shared carry predecessor, not relative to `Full`.")
    lines.append("")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the LUFp benchmark suite and emit a Markdown report.")
    parser.add_argument("--binary", required=True, help="Path to blst_benchmarks")
    parser.add_argument("--profile", choices=sorted(PROFILES.keys()), default="canonical-add-only")
    parser.add_argument("--filter", help="Override benchmark filter regex. Required for exploratory runs.")
    parser.add_argument("--output", help="Markdown report path")
    parser.add_argument("--repetitions", type=int, default=REQUIRED_REPETITIONS)
    parser.add_argument("--min-time", default="0.2s")
    parser.add_argument("--idle-threshold", type=float, default=1.0)
    parser.add_argument("--server-ip", default="", help="Remote benchmark server IP, if any")
    parser.add_argument("--affinity", default="", help="CPU affinity or pinning description")
    args = parser.parse_args()

    if args.profile == "exploratory" and not args.filter:
        raise SystemExit("--filter is required for exploratory profile")

    repo = pathlib.Path(__file__).resolve().parents[1]
    binary = pathlib.Path(args.binary).resolve()
    if not binary.exists():
        raise SystemExit(f"benchmark binary not found: {binary}")

    benchmark_filter = args.filter or PROFILES[args.profile]["filter"]
    if not benchmark_filter:
        raise SystemExit("no benchmark filter available for the selected profile")

    commit, _ = git_metadata(repo)
    stamp = datetime.utcnow().strftime("%Y%m%d-%H%M%S")
    default_output = repo / "docs" / "report" / f"{stamp}-{commit[:12]}-report.md"
    output = pathlib.Path(args.output).resolve() if args.output else default_output

    with tempfile.NamedTemporaryFile(prefix="lufp-bench-", suffix=".json", delete=False) as tmp:
        tmp_path = pathlib.Path(tmp.name)
    try:
        cmd = [
            str(binary),
            f"--benchmark_filter={benchmark_filter}",
            f"--benchmark_repetitions={args.repetitions}",
            f"--benchmark_min_time={args.min_time}",
            "--benchmark_report_aggregates_only=true",
            "--benchmark_out_format=json",
            f"--benchmark_out={tmp_path}",
        ]
        subprocess.run(cmd, cwd=repo, check=True)
        results = parse_benchmark_json(tmp_path)
        render_report(
            repo=repo,
            results=results,
            output_path=output,
            command=" ".join(cmd),
            profile=args.profile,
            repetitions=args.repetitions,
            min_time=args.min_time,
            idle_threshold=args.idle_threshold,
            server_ip=args.server_ip,
            affinity=args.affinity,
        )
        print(output)
    finally:
        tmp_path.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
