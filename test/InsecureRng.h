
struct InsecureRng {
    using u64 = std::uint64_t;

    explicit InsecureRng(u64 seed = 0x9E3779B97F4A7C15ULL) : state(seed ? seed : 1) {}

    // Core generator: xorshift64* (Sebastiano Vigna)
    u64 next_u64() {
        u64 x = state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        state = x;
        return x * 2685821657736338717ULL;
    }

    std::uint64_t uniform_u64() {
		return next_u64();
    }

    // Unbiased integer in [lo, hi] (inclusive), works for lo<=hi.
    std::uint64_t uniform_u64(std::uint64_t lo, std::uint64_t hi) {
        const u64 range = hi - lo + 1;
        const u64 threshold = -range % range;           // rejection to remove modulo bias
        for (;;) {
            u64 r = next_u64();
            if (r >= threshold) return lo + (r % range);
        }
    }

    // Unbiased 32-bit integer in [lo, hi] (inclusive)
    std::uint32_t uniform_u32(std::uint32_t lo, std::uint32_t hi) {
        return static_cast<std::uint32_t>(uniform_u64(lo, hi));
    }

    // Double in [0,1)
    double uniform_double() {
        // Take top 53 bits -> [0, 2^53-1], then scale to [0,1)
        return (next_u64() >> 11) * (1.0 / 9007199254740992.0); // 2^53
    }

private:
    u64 state; // must never be 0; ctor enforces nonzero
};