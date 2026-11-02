#!/usr/bin/env python3

import argparse
import csv
import io
import json
import subprocess


EXECUTABLE = "./out-avx-native/benchmark/middle_product_dispatch_benchmarks"


def run_one(
    mode: str, size: int, algorithm: str, minimum_time: float
) -> dict[str, object]:
    prefix = "one_shot_avx" if mode == "one-shot" else "online_avx"
    name = f"bench_middle_product_{prefix}_{algorithm}_dense/{size}/real_time"
    output = subprocess.check_output(
        [
            "taskset",
            "-c",
            "2",
            EXECUTABLE,
            f"--benchmark_filter=^{name}$",
            f"--benchmark_min_time={minimum_time}s",
            "--benchmark_format=csv",
        ],
        stderr=subprocess.DEVNULL,
        text=True,
    )
    rows = list(csv.DictReader(io.StringIO(output)))
    if len(rows) != 1:
        raise RuntimeError(f"expected one benchmark row for {name}, got {len(rows)}")
    row = rows[0]
    return {
        "size": size,
        "algorithm": algorithm,
        "real_time": float(row["real_time"]),
        "cpu_time": float(row["cpu_time"]),
        "time_unit": row["time_unit"],
        "iterations": int(row["iterations"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("online", "one-shot"), required=True)
    parser.add_argument("--sizes", type=int, nargs="+", required=True)
    parser.add_argument("--repetitions", type=int, default=9)
    parser.add_argument("--min-time", type=float, default=0.2)
    arguments = parser.parse_args()

    samples: list[dict[str, object]] = []
    for repetition in range(arguments.repetitions):
        order = ("fft", "multipair") if repetition % 2 == 0 else ("multipair", "fft")
        for size in arguments.sizes:
            for algorithm in order:
                sample = run_one(
                    arguments.mode, size, algorithm, arguments.min_time
                )
                sample["repetition"] = repetition
                samples.append(sample)
    print(json.dumps(samples, indent=2))


if __name__ == "__main__":
    main()
