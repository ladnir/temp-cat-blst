
#pragma once
#include <vector>



// Standard Bit-Reversal for in-place FFT
template<typename F>
void bit_reverse(std::vector<F>& a) {
    int n = a.size();
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
}

struct DefaultCoeffCtx
{
    template<typename G, typename F>
    static auto plus(const G& a, const F& b) {
        return a + b;
    }
	template<typename G, typename F>
    static auto minus(const G& a, const F& b) {
        return a - b;
    }
    template<typename G, typename F>
    static auto mult(const G& a, const F& b) {
        return a * b;
	}
    template<typename F>
    static auto inverse(const F& a) {
        return a.inverse();
    }

    template<typename F>
    static auto from_u64(std::uint64_t a) {
        return F(a);
    }

};

template <typename G, typename F, typename CoeffCtx = DefaultCoeffCtx>
void fft(std::vector<G>& a, const std::vector<F>& roots, bool invert) {
    int n = a.size();
    bit_reverse(a);

    for (int len = 2; len <= n; len <<= 1) {
        int half = len / 2;
        for (int i = 0; i < n; i += len) {
            for (int j = 0; j < half; j++) {
                // In a real implementation, you'd precompute or select roots 
                // based on 'len' and 'invert'
                F w = roots[j * (n / len)];


                G u = a[i + j];
                G t = CoeffCtx::mult(a[i + j + half], w);
                a[i + j] = CoeffCtx::plus(u, t);
                a[i + j + half] = CoeffCtx::minus(u, t);
            }
        }
    }

    // Integrated scaling for inverse FFT 
    if (invert) {
        F inv_n = CoeffCtx::inverse(CoeffCtx::template from_u64<F>(n));
        for (G& x : a) {
            x = CoeffCtx::mult(x, inv_n); // Scale every group element by 1/n
        }
    }
}