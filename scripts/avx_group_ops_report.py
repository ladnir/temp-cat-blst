#!/usr/bin/env python3
"""Run AVX group-operation benchmarks and render markdown tables.

The runner intentionally executes benchmark passes serially. Do not run a
second copy of this script while one is already benchmarking the same host.
"""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any


CORE_FILTER = (
    r"BM_AvxG1_(Add|AddMixed|Dadd|DaddMixed|Double|Mult)$"
    r"|BM_AvxG2_(Add|AddMixed|Dadd|DaddMixed|Double|Mult)$"
    r"|BM_AvxGT_(Add|Dadd|Double|Mult)$"
    r"|BM_AvxG1_InSubgroup$"
    r"|BM_AvxG2_InSubgroup$"
    r"|BM_AvxG1Affine_InSubgroup$"
    r"|BM_AvxG2Affine_InSubgroup$"
    r"|BM_BlsFp_Sqrt$"
    r"|BM_BlsFp2_Sqrt$"
    r"|BM_AvxPairing_(MillerLoop|FinalExp|Full)$"
    r"|BM_AvxG1_MsmPippenger/(256|512|1024|2048|4096|8192|16384|32768)$"
    r"|BM_AvxG2_MsmPippenger/(256|512|1024|2048|4096|8192|16384|32768)$"
)

BREAKDOWN_FILTER = (
    r"BM_AvxG1_OnCurve$"
    r"|BM_AvxG1_Serialize$"
    r"|BM_AvxG1Affine_OnCurve$"
    r"|BM_AvxG1Affine_Serialize$"
    r"|BM_AvxG1Affine_Deserialize$"
    r"|BM_AvxG2_OnCurve$"
    r"|BM_AvxG2_Serialize$"
    r"|BM_AvxG2Affine_OnCurve$"
    r"|BM_AvxG2Affine_Serialize$"
    r"|BM_AvxG2Affine_Deserialize$"
    r"|BM_AvxGT_InGroup$"
    r"|BM_AvxGT_Serialize$"
    r"|BM_AvxGT_Deserialize$"
    r"|BM_AvxG[12]_Breakdown_ProjectiveToAffine$"
    r"|BM_AvxG[12]_Breakdown_BatchProjectiveToAffine/256$"
    r"|BM_AvxG[12]Affine_Breakdown_(AvxToScalar|ScalarSerialize|ScalarCompress|BytesToScalar|CompressedBytesToScalar|ScalarToAvx|ScalarOnCurve|ScalarInSubgroup)$"
    r"|BM_AvxG[12]Affine_Breakdown_InSubgroupRepeated/(256|1024)$"
    r"|BM_AvxG[12]Affine_Breakdown_BatchInSubgroupMod3Precompute2W40/(256|1024)$"
    r"|BM_AvxGT_Breakdown_(AvxToScalar|ScalarSerialize|BytesToScalar|ScalarToAvx|ScalarInGroup)$"
)


NS_SCALE = {
    "ns": 1.0,
    "us": 1_000.0,
    "ms": 1_000_000.0,
    "s": 1_000_000_000.0,
}


def q_posix(value: str) -> str:
    return shlex.quote(value)


def format_ns(ns: float) -> str:
    if ns >= 1_000:
        return f"{ns / 1_000:,.3f} us"
    if ns >= 100:
        return f"{ns:,.0f} ns"
    return f"{ns:,.1f} ns"


def format_ns_int(ns: float) -> str:
    return f"{ns:,.0f} ns"


def parse_json_files(paths: list[Path]) -> dict[str, dict[str, Any]]:
    rows: dict[str, dict[str, Any]] = {}
    for path in paths:
        with path.open("r", encoding="utf-8-sig") as f:
            data = json.load(f)
        for bench in data.get("benchmarks", []):
            if bench.get("run_type") == "aggregate":
                continue
            unit = bench.get("time_unit", "ns")
            if unit not in NS_SCALE:
                raise ValueError(f"unknown time_unit {unit!r} in {path}")
            bench = dict(bench)
            bench["real_time_ns"] = float(bench["real_time"]) * NS_SCALE[unit]
            bench["cpu_time_ns"] = float(bench["cpu_time"]) * NS_SCALE[unit]
            rows[bench["name"]] = bench
    return rows


def require(rows: dict[str, dict[str, Any]], name: str) -> float:
    try:
        return rows[name]["real_time_ns"]
    except KeyError as exc:
        raise KeyError(f"missing benchmark row: {name}") from exc


def maybe(rows: dict[str, dict[str, Any]], name: str) -> float | None:
    row = rows.get(name)
    if row is None:
        return None
    return row["real_time_ns"]


def run_command(cmd: list[str], stdout_path: Path | None = None) -> None:
    print("+ " + subprocess.list2cmdline(cmd), file=sys.stderr)
    if stdout_path is None:
        subprocess.run(cmd, check=True)
        return
    with stdout_path.open("wb") as f:
        subprocess.run(cmd, check=True, stdout=f)


def build_local_cmd(bench: str, bench_filter: str, out_json: Path, min_time: str) -> list[str]:
    return [
        bench,
        f"--benchmark_filter={bench_filter}",
        f"--benchmark_min_time={min_time}",
        "--benchmark_repetitions=1",
        f"--benchmark_out={str(out_json)}",
        "--benchmark_out_format=json",
    ]


def build_remote_cmd(args: argparse.Namespace, bench_filter: str, out_json: Path, min_time: str) -> list[str]:
    if not args.remote_cwd:
        raise ValueError("--remote-cwd is required with --ssh-host")
    if not args.remote_bench:
        raise ValueError("--remote-bench is required with --ssh-host")

    remote_json = f"/tmp/{out_json.name}"
    remote = (
        f"cd {q_posix(args.remote_cwd)} && "
        f"{q_posix(args.remote_bench)} "
        f"--benchmark_filter={q_posix(bench_filter)} "
        f"--benchmark_min_time={q_posix(min_time)} "
        f"--benchmark_repetitions=1 "
        f"--benchmark_out={q_posix(remote_json)} "
        f"--benchmark_out_format=json >/dev/stderr && "
        f"cat {q_posix(remote_json)}"
    )

    cmd = [
        args.plink,
        "-ssh",
        "-batch",
        "-P",
        str(args.ssh_port),
        "-l",
        args.ssh_user,
    ]
    if args.ssh_key:
        cmd.extend(["-i", args.ssh_key])
    if args.hostkey:
        cmd.extend(["-hostkey", args.hostkey])
    cmd.extend([args.ssh_host, remote])
    return cmd


def run_benchmarks(args: argparse.Namespace, out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    core_json = out_dir / "avx_group_ops_core.json"
    breakdown_json = out_dir / "avx_group_ops_breakdown.json"

    if args.ssh_host:
        run_command(build_remote_cmd(args, CORE_FILTER, core_json, args.min_time), core_json)
        run_command(build_remote_cmd(args, BREAKDOWN_FILTER, breakdown_json, args.min_time), breakdown_json)
    else:
        if not args.bench:
            raise ValueError("--bench is required unless --ssh-host is used")
        run_command(build_local_cmd(args.bench, CORE_FILTER, core_json, args.min_time))
        run_command(build_local_cmd(args.bench, BREAKDOWN_FILTER, breakdown_json, args.min_time))

    return [core_json, breakdown_json]


def emit_core_table(rows: dict[str, dict[str, Any]]) -> list[str]:
    out = ["**Core Ops**", "", "| op | G1 | G2 | GT |", "|---|---:|---:|---:|"]
    entries = [
        ("add `J+J`", "BM_AvxG1_Add", "BM_AvxG2_Add", "BM_AvxGT_Add"),
        ("add `J+A`", "BM_AvxG1_AddMixed", "BM_AvxG2_AddMixed", None),
        ("dadd `J+J`", "BM_AvxG1_Dadd", "BM_AvxG2_Dadd", "BM_AvxGT_Dadd"),
        ("dadd `J+A`", "BM_AvxG1_DaddMixed", "BM_AvxG2_DaddMixed", None),
        ("double", "BM_AvxG1_Double", "BM_AvxG2_Double", "BM_AvxGT_Double"),
        ("mult", "BM_AvxG1_Mult", "BM_AvxG2_Mult", "BM_AvxGT_Mult"),
    ]
    for label, g1, g2, gt in entries:
        cells = [label, format_ns(require(rows, g1)), format_ns(require(rows, g2))]
        cells.append("n/a" if gt is None else format_ns(require(rows, gt)))
        out.append("| " + " | ".join(cells) + " |")
    return out


def emit_validation_table(rows: dict[str, dict[str, Any]]) -> list[str]:
    out = ["", "**Validation / Sqrt**", "", "| op | G1 | G2 | GT / field |", "|---|---:|---:|---:|"]
    entries = [
        ("projective on_curve", "BM_AvxG1_OnCurve", "BM_AvxG2_OnCurve", None),
        ("affine on_curve", "BM_AvxG1Affine_OnCurve", "BM_AvxG2Affine_OnCurve", None),
        ("projective subgroup", "BM_AvxG1_InSubgroup", "BM_AvxG2_InSubgroup", None),
        ("affine subgroup", "BM_AvxG1Affine_InSubgroup", "BM_AvxG2Affine_InSubgroup", None),
    ]
    for label, g1, g2, gt in entries:
        out.append(
            f"| {label} | {format_ns(require(rows, g1))} | {format_ns(require(rows, g2))} | "
            f"{'n/a' if gt is None else format_ns(require(rows, gt))} |"
        )
    out.append(f"| in_group | n/a | n/a | {format_ns(require(rows, 'BM_AvxGT_InGroup'))} |")
    out.append(
        f"| sqrt | `Fp`: {format_ns(require(rows, 'BM_BlsFp_Sqrt'))} | "
        f"`Fp2`: {format_ns(require(rows, 'BM_BlsFp2_Sqrt'))} | n/a |"
    )
    return out


def emit_serde_table(rows: dict[str, dict[str, Any]]) -> list[str]:
    out = ["", "**Serde**", "", "| op | G1 | G2 | GT |", "|---|---:|---:|---:|"]
    entries = [
        ("projective serialize", "BM_AvxG1_Serialize", "BM_AvxG2_Serialize", None),
        ("affine serialize", "BM_AvxG1Affine_Serialize", "BM_AvxG2Affine_Serialize", "BM_AvxGT_Serialize"),
        ("affine deserialize", "BM_AvxG1Affine_Deserialize", "BM_AvxG2Affine_Deserialize", "BM_AvxGT_Deserialize"),
    ]
    for label, g1, g2, gt in entries:
        out.append(
            f"| {label} | {format_ns(require(rows, g1))} | {format_ns(require(rows, g2))} | "
            f"{'n/a' if gt is None else format_ns(require(rows, gt))} |"
        )
    return out


def emit_breakdown_table(rows: dict[str, dict[str, Any]]) -> list[str]:
    def batch_per(name: str) -> float:
        n = int(name.rsplit("/", 1)[1])
        return require(rows, name) / n

    spec = [
        ("G1 projective", "`to_affine`", require(rows, "BM_AvxG1_Breakdown_ProjectiveToAffine")),
        ("G1 projective", "`to_affine_batch_256`", batch_per("BM_AvxG1_Breakdown_BatchProjectiveToAffine/256")),
        ("G1 affine", "AVX -> scalar", require(rows, "BM_AvxG1Affine_Breakdown_AvxToScalar")),
        ("G1 affine", "scalar -> bytes", require(rows, "BM_AvxG1Affine_Breakdown_ScalarSerialize")),
        ("G1 affine", "bytes -> scalar", require(rows, "BM_AvxG1Affine_Breakdown_BytesToScalar")),
        ("G1 affine", "scalar -> AVX", require(rows, "BM_AvxG1Affine_Breakdown_ScalarToAvx")),
        ("G1 affine", "scalar `on_curve`", require(rows, "BM_AvxG1Affine_Breakdown_ScalarOnCurve")),
        ("G1 affine", "scalar `in_subgroup`", require(rows, "BM_AvxG1Affine_Breakdown_ScalarInSubgroup")),
        ("G1 affine", "repeated `in_subgroup`, n=256", batch_per("BM_AvxG1Affine_Breakdown_InSubgroupRepeated/256")),
        ("G1 affine", "batched mod-3 `in_subgroup`, n=256, w=40, precomputed `2P`", batch_per("BM_AvxG1Affine_Breakdown_BatchInSubgroupMod3Precompute2W40/256")),
        ("G1 affine", "repeated `in_subgroup`, n=1024", batch_per("BM_AvxG1Affine_Breakdown_InSubgroupRepeated/1024")),
        ("G1 affine", "batched mod-3 `in_subgroup`, n=1024, w=40, precomputed `2P`", batch_per("BM_AvxG1Affine_Breakdown_BatchInSubgroupMod3Precompute2W40/1024")),
        ("G2 projective", "`to_affine`", require(rows, "BM_AvxG2_Breakdown_ProjectiveToAffine")),
        ("G2 projective", "`to_affine_batch_256`", batch_per("BM_AvxG2_Breakdown_BatchProjectiveToAffine/256")),
        ("G2 affine", "AVX -> scalar", require(rows, "BM_AvxG2Affine_Breakdown_AvxToScalar")),
        ("G2 affine", "scalar -> bytes", require(rows, "BM_AvxG2Affine_Breakdown_ScalarSerialize")),
        ("G2 affine", "bytes -> scalar", require(rows, "BM_AvxG2Affine_Breakdown_BytesToScalar")),
        ("G2 affine", "scalar -> AVX", require(rows, "BM_AvxG2Affine_Breakdown_ScalarToAvx")),
        ("G2 affine", "scalar `on_curve`", require(rows, "BM_AvxG2Affine_Breakdown_ScalarOnCurve")),
        ("G2 affine", "scalar `in_subgroup`", require(rows, "BM_AvxG2Affine_Breakdown_ScalarInSubgroup")),
        ("G2 affine", "repeated `in_subgroup`, n=256", batch_per("BM_AvxG2Affine_Breakdown_InSubgroupRepeated/256")),
        ("G2 affine", "batched mod-3 `in_subgroup`, n=256, w=40, precomputed `2P`", batch_per("BM_AvxG2Affine_Breakdown_BatchInSubgroupMod3Precompute2W40/256")),
        ("G2 affine", "repeated `in_subgroup`, n=1024", batch_per("BM_AvxG2Affine_Breakdown_InSubgroupRepeated/1024")),
        ("G2 affine", "batched mod-3 `in_subgroup`, n=1024, w=40, precomputed `2P`", batch_per("BM_AvxG2Affine_Breakdown_BatchInSubgroupMod3Precompute2W40/1024")),
        ("GT", "AVX -> scalar", require(rows, "BM_AvxGT_Breakdown_AvxToScalar")),
        ("GT", "scalar -> bytes", require(rows, "BM_AvxGT_Breakdown_ScalarSerialize")),
        ("GT", "bytes -> scalar", require(rows, "BM_AvxGT_Breakdown_BytesToScalar")),
        ("GT", "scalar -> AVX", require(rows, "BM_AvxGT_Breakdown_ScalarToAvx")),
        ("GT", "scalar `in_group`", require(rows, "BM_AvxGT_Breakdown_ScalarInGroup")),
    ]

    out = ["", "**Serialize / Deserialize Micro Breakdown**", "", "| Group | Operation / Step | Per Point |", "|---|---|---:|"]
    for group, op, ns in spec:
        out.append(f"| {group} | {op} | {format_ns_int(ns)} |")
    return out


def emit_projective_path_table(rows: dict[str, dict[str, Any]]) -> list[str]:
    g1_single = (
        require(rows, "BM_AvxG1_Breakdown_ProjectiveToAffine")
        + require(rows, "BM_AvxG1Affine_Breakdown_AvxToScalar")
        + require(rows, "BM_AvxG1Affine_Breakdown_ScalarSerialize")
    )
    g1_batch = (
        require(rows, "BM_AvxG1_Breakdown_BatchProjectiveToAffine/256") / 256
        + require(rows, "BM_AvxG1Affine_Breakdown_AvxToScalar")
        + require(rows, "BM_AvxG1Affine_Breakdown_ScalarSerialize")
    )
    g2_single = (
        require(rows, "BM_AvxG2_Breakdown_ProjectiveToAffine")
        + require(rows, "BM_AvxG2Affine_Breakdown_AvxToScalar")
        + require(rows, "BM_AvxG2Affine_Breakdown_ScalarSerialize")
    )
    g2_batch = (
        require(rows, "BM_AvxG2_Breakdown_BatchProjectiveToAffine/256") / 256
        + require(rows, "BM_AvxG2Affine_Breakdown_AvxToScalar")
        + require(rows, "BM_AvxG2Affine_Breakdown_ScalarSerialize")
    )

    out = ["", "**Projective Serialize Path, Per Point**", "", "| Group | Path | Per Point |", "|---|---|---:|"]
    out.append(f"| G1 | single `to_affine` + AVX -> scalar + scalar -> bytes | {format_ns_int(g1_single)} |")
    out.append(f"| G1 | batch-256 `to_affine` + AVX -> scalar + scalar -> bytes | {format_ns_int(g1_batch)} |")
    out.append(f"| G2 | single `to_affine` + AVX -> scalar + scalar -> bytes | {format_ns_int(g2_single)} |")
    out.append(f"| G2 | batch-256 `to_affine` + AVX -> scalar + scalar -> bytes | {format_ns_int(g2_batch)} |")
    return out


def emit_msm_table(rows: dict[str, dict[str, Any]]) -> list[str]:
    out = ["", "**MSM, Amortized Per Term**", "", "| n | G1 | G2 |", "|---:|---:|---:|"]
    for n in [256, 512, 1024, 2048, 4096, 8192, 16384, 32768]:
        g1 = require(rows, f"BM_AvxG1_MsmPippenger/{n}") / n
        g2 = require(rows, f"BM_AvxG2_MsmPippenger/{n}") / n
        out.append(f"| {n} | {format_ns(g1)} | {format_ns(g2)} |")
    return out


def emit_pairing_table(rows: dict[str, dict[str, Any]]) -> list[str]:
    out = ["", "**Pairing**", "", "| op | time |", "|---|---:|"]
    out.append(f"| Miller loop | {format_ns(require(rows, 'BM_AvxPairing_MillerLoop'))} |")
    out.append(f"| Final exponentiation | {format_ns(require(rows, 'BM_AvxPairing_FinalExp'))} |")
    out.append(f"| Full pairing | {format_ns(require(rows, 'BM_AvxPairing_Full'))} |")
    return out


def emit_report(rows: dict[str, dict[str, Any]]) -> str:
    lines: list[str] = []
    lines.extend(emit_core_table(rows))
    lines.extend(emit_validation_table(rows))
    lines.extend(emit_serde_table(rows))
    lines.extend(emit_breakdown_table(rows))
    lines.extend(emit_projective_path_table(rows))
    lines.extend(emit_msm_table(rows))
    lines.extend(emit_pairing_table(rows))
    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bench", help="local path to avx_group_ops_benchmarks")
    parser.add_argument("--out-dir", default="out/avx_group_ops_report", help="directory for JSON and markdown output")
    parser.add_argument("--report", help="markdown output path; defaults to <out-dir>/report.md")
    parser.add_argument("--min-time", default="0.25s", help="google-benchmark --benchmark_min_time value")
    parser.add_argument("--from-json", nargs="+", type=Path, help="render from existing google-benchmark JSON files")

    parser.add_argument("--ssh-host", help="run benchmark remotely through plink")
    parser.add_argument("--ssh-port", type=int, default=22)
    parser.add_argument("--ssh-user")
    parser.add_argument("--ssh-key")
    parser.add_argument("--hostkey", help="plink -hostkey value")
    parser.add_argument("--plink", default="plink")
    parser.add_argument("--remote-cwd", help="remote working directory")
    parser.add_argument("--remote-bench", help="remote path to avx_group_ops_benchmarks")

    args = parser.parse_args()
    if args.ssh_host and not args.ssh_user:
        parser.error("--ssh-user is required with --ssh-host")
    return args


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    report_path = Path(args.report) if args.report else out_dir / "report.md"

    if args.from_json:
        json_paths = args.from_json
    else:
        json_paths = run_benchmarks(args, out_dir)

    rows = parse_json_files(json_paths)
    report = emit_report(rows)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(report, encoding="utf-8")
    print(report)
    print(f"Wrote {report_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
