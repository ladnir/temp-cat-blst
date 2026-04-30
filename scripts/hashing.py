import argparse
import math


PAPER_BATCHES = (32, 64, 128, 256, 512)
PAPER_D_VALUES = tuple(range(8, 17))
PAPER_M_VALUES = (8, 16, 32, 64, 128, 256, 512, 1024)


CHOSEN_40_BIT_C2 = {
    32: {"copies": 2, "d": 13, "m": 8, "supported": 38},
    64: {"copies": 2, "d": 16, "m": 8, "supported": 69},
    128: {"copies": 2, "d": 14, "m": 16, "supported": 162},
    256: {"copies": 2, "d": 12, "m": 32, "supported": 351},
    512: {"copies": 2, "d": 10, "m": 64, "supported": 606},
}

CHOSEN_128_BIT_C20 = {
    32: {"copies": 20, "d": 8, "m": 8, "bits": 512.3},
    64: {"copies": 20, "d": 9, "m": 8, "bits": 243.6},
    128: {"copies": 20, "d": 9, "m": 16, "bits": 475.5},
    256: {"copies": 20, "d": 9, "m": 32, "bits": 936.6},
    512: {"copies": 20, "d": 9, "m": 64, "bits": 1125.9},
}


def _dw_fj(j, m, d, epsilon):
    item_count = m * d / (1 + epsilon)
    if j <= 0 or j > m:
        return 0.0

    log_binom = math.lgamma(m + 1) - math.lgamma(j + 1) - math.lgamma(m - j + 1)
    ratio1 = item_count * ((j / m) ** 2) / (j * d)
    if ratio1 <= 0.0:
        return 0.0
    log_term2 = (j * d) * math.log(ratio1)

    if item_count - j * d == 0:
        log_term3 = 0.0
    else:
        ratio2 = (item_count - item_count * ((j / m) ** 2)) / (item_count - j * d)
        if ratio2 <= 0.0:
            return float("inf")
        log_term3 = (item_count - j * d) * math.log(ratio2)

    return math.exp(log_binom + log_term2 + log_term3)


def dw_c2_failure_bound(m, d, item_count):
    epsilon = (m * d / item_count) - 1
    upper = math.floor(m / (1 + epsilon))
    return sum(_dw_fj(j, m, d, epsilon) for j in range(1, upper + 1))


def max_items_dw_c2(m, d, target_bits):
    fail_prob = 2 ** (-target_bits)
    lo = 1
    hi = m * d
    best = 0
    while lo <= hi:
        mid = (lo + hi) // 2
        if dw_c2_failure_bound(m, d, mid) < fail_prob:
            best = mid
            lo = mid + 1
        else:
            hi = mid - 1
    return best


def _logsumexp(values):
    values = [value for value in values if value != -math.inf]
    if not values:
        return -math.inf
    max_value = max(values)
    return max_value + math.log(sum(math.exp(value - max_value) for value in values))


def _log_binom(n, k):
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _log_binomial_tail(n, p, lower):
    if lower > n:
        return -math.inf
    if p <= 0.0:
        return -math.inf
    if p >= 1.0:
        return 0.0

    log_p = math.log(p)
    log_q = math.log1p(-p)
    return _logsumexp(
        _log_binom(n, k) + k * log_p + (n - k) * log_q
        for k in range(lower, n + 1)
    )


def heuristic_failure_bits(item_count, m, d, copies):
    """Union-bound Hall obstructions for arbitrary copy count.

    For each bucket set S of size j, a transaction is trapped inside S if all
    copies/hash choices land in S. If more than d*j transactions are trapped,
    those buckets cannot be filled. This is the heuristic 128-bit criterion used
    for the paper's high-copy BEAT++ comparison.
    """
    log_terms = []
    max_j = min(m, (item_count - 1) // d)
    for j in range(1, max_j + 1):
        trapped_probability = (j / m) ** copies
        log_terms.append(
            _log_binom(m, j)
            + _log_binomial_tail(item_count, trapped_probability, d * j + 1)
        )

    log_fail = _logsumexp(log_terms)
    if log_fail == -math.inf:
        return float("inf")
    return -log_fail / math.log(2.0)


def best_c2_parameters(batch_size, target_bits):
    candidates = []
    for d in PAPER_D_VALUES:
        for m in PAPER_M_VALUES:
            if d * m <= batch_size:
                continue
            supported = max_items_dw_c2(m, d, target_bits)
            if supported > batch_size:
                candidates.append((d * m, d * m * math.log2(m), d, m, supported))
    if not candidates:
        raise RuntimeError(f"no c=2 parameters support batch {batch_size}")
    slots, _, d, m, supported = min(candidates)
    return {"copies": 2, "d": d, "m": m, "supported": supported, "slots": slots}


def best_heuristic_parameters(batch_size, copies, target_bits):
    candidates = []
    for d in PAPER_D_VALUES:
        for m in PAPER_M_VALUES:
            if d * m <= batch_size:
                continue
            bits = heuristic_failure_bits(batch_size, m, d, copies)
            if bits >= target_bits:
                candidates.append((d * m, d * m * math.log2(m), d, m, bits))
    if not candidates:
        raise RuntimeError(f"no c={copies} heuristic parameters support batch {batch_size}")
    slots, _, d, m, bits = min(candidates)
    return {"copies": copies, "d": d, "m": m, "bits": bits, "slots": slots}


def print_table(title, rows):
    print(title)
    print("B,c,d,m,slots,supported_or_bits,load")
    for batch_size in PAPER_BATCHES:
        row = rows[batch_size]
        slots = row["d"] * row["m"]
        support = row.get("supported", row.get("bits"))
        if isinstance(support, float):
            support = f"{support:.1f}"
        print(
            f"{batch_size},{row['copies']},{row['d']},{row['m']},"
            f"{slots},{support},{batch_size / slots:.3f}"
        )
    print()


def main():
    parser = argparse.ArgumentParser(description="Search BEAT++ cuckoo-hashing parameters.")
    parser.add_argument("--target-bits", type=int, default=40)
    parser.add_argument("--copies", type=int, default=20)
    args = parser.parse_args()

    c2_rows = {
        batch_size: best_c2_parameters(batch_size, args.target_bits)
        for batch_size in PAPER_BATCHES
    }
    heuristic_rows = {
        batch_size: best_heuristic_parameters(batch_size, args.copies, 128)
        for batch_size in PAPER_BATCHES
    }

    print_table(f"best c=2, {args.target_bits}-bit DW-style bound", c2_rows)
    print_table(f"best c={args.copies}, 128-bit Hall-obstruction heuristic", heuristic_rows)
    print_table("chosen paper c=2, 40-bit set", CHOSEN_40_BIT_C2)
    print_table("chosen paper c=20, 128-bit heuristic set", CHOSEN_128_BIT_C20)


if __name__ == "__main__":
    main()
