#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import os
import pathlib
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Tuple


SUITE_FILTER = (
    r"^(BM_UFpAddChain.*|"
    r"BM_UFpAddAccounting.*|"
    r"BM_Fp6MulUfpWide.*|"
    r"BM_Fp6MulUfpAccounting.*|"
    r"BM_LUFpAddChain.*|"
    r"BM_LUFpAddAccounting.*|"
    r"BM_Fp6MulLufpWide.*|"
    r"BM_Fp6MulLufpAccounting.*|"
    r"BM_LUFpNormalizeHiNeg1Pos2.*|"
    r"BM_LUFpCarryChain.*|"
    r"BM_LUFpRepairPrefixLex.*)$"
)

EXPECTED_FAMILIES = {
    "UFp add": ["BM_UFpAddChain/256"],
    "UFp add accounting": [
        "BM_UFpAddAccounting<UFpAddStage::RawAdd>/256",
        "BM_UFpAddAccounting<UFpAddStage::RawAddBias>/256",
        "BM_UFpAddAccounting<UFpAddStage::RawAddBiasCarry>/256",
        "BM_UFpAddAccounting<UFpAddStage::Full>/256",
    ],
    "UFp Fp6 mul": ["BM_Fp6MulUfpWide/256"],
    "UFp Fp6 mul accounting": [
        "BM_Fp6MulUfpAccounting<UFpFp6MulStage::PreAdd>/256",
        "BM_Fp6MulUfpAccounting<UFpFp6MulStage::PreAddMul>/256",
        "BM_Fp6MulUfpAccounting<UFpFp6MulStage::Full>/256",
    ],
    "LUFp add": ["BM_LUFpAddChain/256"],
    "LUFp add accounting": [
        "BM_LUFpAddAccounting<LUFpAddStage::RawAdd>/256",
        "BM_LUFpAddAccounting<LUFpAddStage::RawAddBias>/256",
        "BM_LUFpAddAccounting<LUFpAddStage::RawAddBiasCarry>/256",
        "BM_LUFpAddAccounting<LUFpAddStage::Full>/256",
    ],
    "LUFp renorm": [
        "BM_LUFpNormalizeHiNeg1Pos2/256",
        "BM_LUFpCarryChain/256",
        "BM_LUFpRepairPrefixLex<1>/256",
        "BM_LUFpRepairPrefixLex<2>/256",
        "BM_LUFpRepairPrefixLex<4>/256",
    ],
    "LUFp Fp6 mul": ["BM_Fp6MulLufpWide/256"],
    "LUFp Fp6 mul accounting": [
        "BM_Fp6MulLufpAccounting<LUFpFp6MulStage::PreAdd>/256",
        "BM_Fp6MulLufpAccounting<LUFpFp6MulStage::PreAddMul>/256",
        "BM_Fp6MulLufpAccounting<LUFpFp6MulStage::PreAddMulRawCombine>/256",
        "BM_Fp6MulLufpAccounting<LUFpFp6MulStage::Full>/256",
    ],
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
        if not agg:
            continue
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


def collect_red_flags(
    results: Dict[str, Aggregate],
    idle_threshold: float,
    load1: float,
) -> List[Tuple[str, str, str, str]]:
    flags: List[Tuple[str, str, str, str]] = []

    for family, names in EXPECTED_FAMILIES.items():
        missing = [name for name in names if name not in results]
        if missing:
            flags.append(
                (
                    "fatal",
                    family,
                    f"missing benchmarks: {', '.join(missing)}",
                    "fix benchmark registration or filter before trusting the run",
                )
            )

    for name, agg in results.items():
        if agg.cv > 5.0:
            flags.append(("fatal", name, f"coefficient of variation is {fmt_pct(agg.cv)}", "rerun on a more stable idle machine"))
        elif agg.cv > 2.0:
            flags.append(("warning", name, f"coefficient of variation is {fmt_pct(agg.cv)}", "rerun to improve stability"))

    if load1 > idle_threshold * 2.0:
        flags.append(("fatal", "environment", f"load average before run was {load1:.2f} (threshold {idle_threshold:.2f})", "rerun on an idle machine"))
    elif load1 > idle_threshold:
        flags.append(("warning", "environment", f"load average before run was {load1:.2f} (threshold {idle_threshold:.2f})", "prefer an idler machine or rerun"))

    return flags


def accounting_rows(results: Dict[str, Aggregate], names: List[str]) -> List[Tuple[str, float, float, float]]:
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


def op_summary(results: Dict[str, Aggregate]) -> List[Tuple[str, str, str]]:
    items: List[Tuple[str, str, str]] = []

    ufp_add = results["BM_UFpAddAccounting<UFpAddStage::Full>/256"].mean
    lufp_add = results["BM_LUFpAddAccounting<LUFpAddStage::Full>/256"].mean
    add_ops = 256 * 16
    items.append(("UFp add", "1 UFp add", fmt_ns(ufp_add / add_ops)))
    items.append(("LUFp add", "1 LUFp add + exact renorm", fmt_ns(lufp_add / add_ops)))

    ufp_fp6 = accounting_rows(
        results,
        [
            "BM_Fp6MulUfpAccounting<UFpFp6MulStage::PreAdd>/256",
            "BM_Fp6MulUfpAccounting<UFpFp6MulStage::PreAddMul>/256",
            "BM_Fp6MulUfpAccounting<UFpFp6MulStage::Full>/256",
        ],
    )
    items.append(("UFp Fp6 mul", "narrow pre-add (6 Fp2 adds)", fmt_ns(ufp_fp6[0][2] / (256 * 6))))
    items.append(("UFp Fp6 mul", "Fp2 mul (6 muls)", fmt_ns(ufp_fp6[1][2] / (256 * 6))))
    items.append(("UFp Fp6 mul", "wide combine add/sub (26 UFp ops)", fmt_ns(ufp_fp6[2][2] / (256 * 26))))

    lufp_fp6 = accounting_rows(
        results,
        [
            "BM_Fp6MulLufpAccounting<LUFpFp6MulStage::PreAdd>/256",
            "BM_Fp6MulLufpAccounting<LUFpFp6MulStage::PreAddMul>/256",
            "BM_Fp6MulLufpAccounting<LUFpFp6MulStage::PreAddMulRawCombine>/256",
            "BM_Fp6MulLufpAccounting<LUFpFp6MulStage::Full>/256",
        ],
    )
    items.append(("LUFp Fp6 mul", "narrow pre-add (6 Fp2 adds)", fmt_ns(lufp_fp6[0][2] / (256 * 6))))
    items.append(("LUFp Fp6 mul", "Fp2 mul (6 muls)", fmt_ns(lufp_fp6[1][2] / (256 * 6))))
    items.append(("LUFp Fp6 mul", "raw lazy combine (26 UFp ops)", fmt_ns(lufp_fp6[2][2] / (256 * 26))))
    items.append(("LUFp Fp6 mul", "renorm/output boundary (6 outputs)", fmt_ns(lufp_fp6[3][2] / (256 * 6))))

    return items


def speed_of_light_section(results: Dict[str, Aggregate]) -> List[str]:
    lines: List[str] = []
    families = [
        (
            "UFp add",
            "BM_UFpAddAccounting<UFpAddStage::RawAdd>/256",
            "BM_UFpAddAccounting<UFpAddStage::Full>/256",
            "Lower bound is the measured raw-add stage. Assumption: the full kernel must include that stage and all later stages are non-negative overhead.",
        ),
        (
            "LUFp add",
            "BM_LUFpAddAccounting<LUFpAddStage::RawAdd>/256",
            "BM_LUFpAddAccounting<LUFpAddStage::Full>/256",
            "Lower bound is the measured raw lazy add stage. Assumption: bias, carry, and repair are non-negative overhead beyond raw add.",
        ),
        (
            "UFp Fp6 mul",
            "BM_Fp6MulUfpAccounting<UFpFp6MulStage::PreAddMul>/256",
            "BM_Fp6MulUfpAccounting<UFpFp6MulStage::Full>/256",
            "Lower bound is the cumulative pre-add plus multiply stage. Assumption: wide combine and output handling are non-negative overhead.",
        ),
        (
            "LUFp Fp6 mul",
            "BM_Fp6MulLufpAccounting<LUFpFp6MulStage::PreAddMulRawCombine>/256",
            "BM_Fp6MulLufpAccounting<LUFpFp6MulStage::Full>/256",
            "Lower bound is the cumulative pre-add plus multiply plus raw-lazy-combine stage. Assumption: the final boundary is non-negative overhead.",
        ),
    ]

    for label, lower_name, measured_name, assumptions in families:
        lower = results[lower_name].mean
        measured = results[measured_name].mean
        ratio = 0.0 if lower == 0.0 else measured / lower
        gap = measured - lower
        lines.append(f"### {label}")
        lines.append(f"- Lower bound: {fmt_ns(lower)}")
        lines.append(f"- Measured: {fmt_ns(measured)}")
        lines.append(f"- Gap: {fmt_ns(gap)} ({ratio:.2f}x of lower bound)")
        lines.append(f"- Assumptions: {assumptions}")
        lines.append(f"- Facts: lower and measured values come directly from the accounting suite for this run.")
        lines.append(f"- Inferences: the gap is attributable to stages not included in the chosen lower-bound stage.")
        lines.append(f"- Speculation: none beyond the measured stage decomposition.")
        lines.append("")
    return lines


def render_report(
    repo: pathlib.Path,
    results: Dict[str, Aggregate],
    command: str,
    output_path: pathlib.Path,
    idle_threshold: float,
) -> None:
    timestamp = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")
    commit, dirty = git_metadata(repo)
    load1, load5, load15 = load_average()
    flags = collect_red_flags(results, idle_threshold, load1)

    # Integrity check: accounting full should match the top-level full benchmark closely.
    comparisons = [
        ("UFp add", "BM_UFpAddChain/256", "BM_UFpAddAccounting<UFpAddStage::Full>/256"),
        ("UFp Fp6 mul", "BM_Fp6MulUfpWide/256", "BM_Fp6MulUfpAccounting<UFpFp6MulStage::Full>/256"),
        ("LUFp Fp6 mul", "BM_Fp6MulLufpWide/256", "BM_Fp6MulLufpAccounting<LUFpFp6MulStage::Full>/256"),
    ]
    for family, full_name, accounting_name in comparisons:
        if full_name in results and accounting_name in results:
            full = results[full_name].mean
            accounting = results[accounting_name].mean
            if full != 0.0 and abs(accounting - full) / full > 3.0 / 100.0:
                flags.append(
                    (
                        "fatal",
                        family,
                        f"accounting full {fmt_ns(accounting)} does not match full benchmark {fmt_ns(full)} within 3%",
                        "align the full benchmark and accounting Full stage before trusting the decomposition",
                    )
                )

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
    lines.append(f"- Git commit: `{commit}`")
    lines.append(f"- Git status: `{dirty}`")
    lines.append(f"- Hostname: `{os.uname().nodename}`")
    lines.append(f"- CPU model: {cpu_model()}")
    lines.append(f"- Load average before run: {load1:.2f}, {load5:.2f}, {load15:.2f}")
    lines.append(f"- Idle threshold: {idle_threshold:.2f}")
    lines.append(f"- Benchmark command: `{command}`")
    lines.append("")

    lines.append("## Benchmark Coverage")
    for family, names in EXPECTED_FAMILIES.items():
        missing = [name for name in names if name not in results]
        status = "present" if not missing else "missing"
        lines.append(f"- {family}: {status}")
        if missing:
            lines.append(f"  missing: {', '.join(missing)}")
    lines.append("")

    lines.append("## Raw Results")
    raw_rows = []
    for name in sorted(results.keys()):
        agg = results[name]
        raw_rows.append([name, fmt_ns(agg.mean), fmt_ns(agg.median), fmt_ns(agg.stddev), fmt_pct(agg.cv)])
    add_table(lines, ["Benchmark", "Mean", "Median", "Stddev", "CV"], raw_rows)

    lines.append("## Speed-of-Light Summary")
    lines.extend(speed_of_light_section(results))

    lines.append("## Accounting Breakdown")
    accounting_families = {
        "UFp add": EXPECTED_FAMILIES["UFp add accounting"],
        "UFp Fp6 mul": EXPECTED_FAMILIES["UFp Fp6 mul accounting"],
        "LUFp add": EXPECTED_FAMILIES["LUFp add accounting"],
        "LUFp Fp6 mul": EXPECTED_FAMILIES["LUFp Fp6 mul accounting"],
    }
    for label, names in accounting_families.items():
        if any(name not in results for name in names):
            continue
        lines.append(f"### {label}")
        rows = [
            [name, fmt_ns(cumulative), fmt_ns(delta), fmt_pct(pct)]
            for name, cumulative, delta, pct in accounting_rows(results, names)
        ]
        add_table(lines, ["Stage", "Cumulative", "Delta", "% Full"], rows)

    lines.append("## Per-Operation Summary")
    per_op_rows = [[family, op, value] for family, op, value in op_summary(results)]
    add_table(lines, ["Family", "Operation", "Average"], per_op_rows)

    lines.append("## Full Comparisons")
    full_rows = []
    if "BM_UFpAddChain/256" in results and "BM_LUFpAddAccounting<LUFpAddStage::Full>/256" in results:
        a = results["BM_UFpAddChain/256"].mean
        b = results["BM_LUFpAddAccounting<LUFpAddStage::Full>/256"].mean
        full_rows.append(["UFp add vs LUFp add", fmt_ns(a), fmt_ns(b), fmt_ns(b - a), fmt_pct(((b - a) / a) * 100.0)])
    if "BM_Fp6MulUfpWide/256" in results and "BM_Fp6MulLufpWide/256" in results:
        a = results["BM_Fp6MulUfpWide/256"].mean
        b = results["BM_Fp6MulLufpWide/256"].mean
        full_rows.append(["UFp Fp6 mul vs LUFp Fp6 mul", fmt_ns(a), fmt_ns(b), fmt_ns(b - a), fmt_pct(((b - a) / a) * 100.0)])
    add_table(lines, ["Comparison", "Baseline", "Candidate", "Absolute Delta", "Percent Delta"], full_rows)

    lines.append("## Interpretation")
    lines.append("- Facts: this report uses the accounting benchmark families as the source of constituent timing, and the top-level benchmark families as the direct end-to-end comparison points.")
    lines.append("- Inferences: any stage gap reported above is only attributable to stages excluded by the chosen lower-bound stage.")
    lines.append("- Speculation: none beyond the explicitly labeled speed-of-light sections.")
    lines.append("")

    lines.append("## Notes")
    lines.append("- The speed-of-light numbers in this report are measured structural lower bounds from the accounting suite, not hardware roofline estimates.")
    lines.append("- Any fatal red flag at the top means the run should not be used for conclusions until the flagged issue is fixed.")
    lines.append("")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the LUFp benchmark suite and emit a Markdown report.")
    parser.add_argument("--binary", required=True, help="Path to blst_benchmarks")
    parser.add_argument("--output", help="Markdown report path")
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--min-time", default="0.05s")
    parser.add_argument("--idle-threshold", type=float, default=1.0)
    args = parser.parse_args()

    repo = pathlib.Path(__file__).resolve().parents[1]
    binary = pathlib.Path(args.binary).resolve()
    if not binary.exists():
        raise SystemExit(f"benchmark binary not found: {binary}")

    commit, _ = git_metadata(repo)
    stamp = datetime.utcnow().strftime("%Y%m%d-%H%M%S")
    default_output = repo / "docs" / "report" / f"{stamp}-{commit[:12]}-report.md"
    output = pathlib.Path(args.output).resolve() if args.output else default_output

    with tempfile.NamedTemporaryFile(prefix="lufp-bench-", suffix=".json", delete=False) as tmp:
        tmp_path = pathlib.Path(tmp.name)
    try:
        cmd = [
            str(binary),
            f"--benchmark_filter={SUITE_FILTER}",
            f"--benchmark_repetitions={args.repetitions}",
            f"--benchmark_min_time={args.min_time}",
            "--benchmark_report_aggregates_only=true",
            "--benchmark_out_format=json",
            f"--benchmark_out={tmp_path}",
        ]
        subprocess.run(cmd, cwd=repo, check=True)
        results = parse_benchmark_json(tmp_path)
        render_report(repo, results, " ".join(cmd), output, args.idle_threshold)
        print(output)
    finally:
        tmp_path.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
