#!/usr/bin/env python3

import json
import random
import statistics
import sys


def percentile(sorted_values: list[float], probability: float) -> float:
    return sorted_values[round(probability * (len(sorted_values) - 1))]


def main() -> None:
    with open(sys.argv[1], encoding="utf-8") as source:
        samples = json.load(source)

    rng = random.Random(0xC20550)
    for size in sorted({sample["size"] for sample in samples}):
        by_repetition: dict[int, dict[str, float]] = {}
        for sample in samples:
            if sample["size"] == size:
                by_repetition.setdefault(sample["repetition"], {})[
                    sample["algorithm"]
                ] = sample["real_time"]

        fft = [values["fft"] for values in by_repetition.values()]
        multipair = [values["multipair"] for values in by_repetition.values()]
        paired_percent = [
            100.0 * (multipair_time - fft_time) / fft_time
            for fft_time, multipair_time in zip(fft, multipair)
        ]
        bootstrap_medians = sorted(
            statistics.median(rng.choices(paired_percent, k=len(paired_percent)))
            for _ in range(100_000)
        )
        print(
            f"B={size}: fft={statistics.median(fft):.6f} ms, "
            f"multipair={statistics.median(multipair):.6f} ms, "
            f"paired_delta={statistics.median(paired_percent):+.3f}% "
            f"(bootstrap95 {percentile(bootstrap_medians, 0.025):+.3f}%.."
            f"{percentile(bootstrap_medians, 0.975):+.3f}%)"
        )


if __name__ == "__main__":
    main()
