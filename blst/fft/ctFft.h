// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
//
// CatBLST radix-two FFT and shape-planned truncated-transform primitives.

#pragma once
#include <algorithm>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string_view>
#include <vector>
#include "blst/common/defines.h"
#include "blst/common/timer.h"

namespace blst::fft_detail
{
    [[nodiscard]] constexpr bool valid_transform_size(
        const std::size_t n) noexcept
    {
        return n != 0 && (n & (n - 1)) == 0 &&
            n <= static_cast<std::size_t>(std::numeric_limits<int>::max());
    }
}

// Applies the binary bit-reversal permutation in place. The span length must
// fit the signed-index transform implementation; FFT callers additionally
// require a nonzero power-of-two length.
template<typename F>
void bit_reverse(std::span<F> a) {
	if (a.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
		throw std::invalid_argument("FFT size exceeds the supported index range");
	}
	const int n = static_cast<int>(a.size());
	for (int i = 1, j = 0; i < n; i++) {
		int bit = n >> 1;
		for (; j & bit; bit >>= 1) j ^= bit;
		j ^= bit;
		if (i < j) std::swap(a[i], a[j]);
	}
}

// In-place radix-two DIT butterflies. `a.size()` must be a supported nonzero
// power of two. `roots` must contain at least that many consecutive powers of
// the corresponding primitive root, starting with one, and must not overlap
// `a`. The schedule and root accesses depend only on the public transform size.
template <bool InverseFromForwardRoots, typename G, typename F, typename CoeffCtx>
void fft_butterflies(std::span<G> a, std::span<const F> roots) {
	BLST_ASSERT(blst::fft_detail::valid_transform_size(a.size()));
	BLST_ASSERT(roots.size() >= a.size());
	bit_reverse(a);
	const int n = static_cast<int>(a.size());

	for (int len = 2; len <= n;) {
		int half = len / 2;
		for (int i = 0; i < n; i += len) {

			// j == 0 case, this is 20% of the work
			// and we get to skip the multiplication
			// because w == 1.
			{ 
				G u = a[i];
				G t = a[i + half];
				a[i] = CoeffCtx::plus(u, t);
				a[i + half] = CoeffCtx::minus(u, t);
			}

			for (int j = 1; j < half; j++) {
				int root_index = j * (n / len);
				if constexpr (InverseFromForwardRoots) {
					root_index = n - root_index;
				}
				const F& w = roots[root_index];

				G u = a[i + j];
				G t = CoeffCtx::mult(a[i + j + half], w);
				a[i + j] = CoeffCtx::plus(u, t);
				a[i + j + half] = CoeffCtx::minus(u, t);
			}
		}
		if (len == n) break;
		len <<= 1;
	}
}

template <typename G, typename F, typename CoeffCtx>
void fft_scale(std::span<G> a, const F& inverse_size) {
	for (G& x : a) {
		x = CoeffCtx::mult(x, inverse_size);
	}
}

// Runs an in-place transform using a direction-specific root table. Forward
// calls require forward powers. Inverse calls require inverse powers and apply
// the 1/n scale. `invert` does not change the supplied root direction.
template <typename G, typename F, typename CoeffCtx>
void fft(std::span<G> a, std::span<const F> roots, bool invert) {
	fft_butterflies<false, G, F, CoeffCtx>(a, roots);

	// Integrated scaling for inverse FFT.
	if (invert) {
		F inv_n = CoeffCtx::inverse(CoeffCtx::template from_u64<F>(a.size()));
		fft_scale<G, F, CoeffCtx>(a, inv_n);
	}
}

// Runs an in-place transform from one forward-root table. The inverse path
// reads reciprocal powers by index and applies the caller-prepared 1/n scale.
// The data, root table, and scale object must not overlap.
template <typename G, typename F, typename CoeffCtx>
void fft_with_forward_roots(
	std::span<G> a,
	std::span<const F> forward_roots,
	bool invert,
	const F& inverse_size) {
	if (invert) {
		fft_butterflies<true, G, F, CoeffCtx>(a, forward_roots);
		fft_scale<G, F, CoeffCtx>(a, inverse_size);
	} else {
		fft_butterflies<false, G, F, CoeffCtx>(a, forward_roots);
	}
}

// Computes the complete forward transform when only `a[0..xn)` may be
// nonidentity. Every element in `a[xn..n)` must be the additive identity. The
// plan depends only on public sizes, may be shared by concurrent readers, and
// must have been created for exactly the supplied `(n, xn)` shape. Data and
// roots must not overlap. All n outputs are valid on return.
template <typename G, typename F, typename CoeffCtx>
void fft_forward_truncated_prefix(
        std::span<G> a,
        std::span<const F> roots,
        std::size_t xn,
        const struct FftTruncatedForwardPlan& plan,
        blst::TimingStats* timings = nullptr,
        std::string_view metric_prefix = "fft_trunc");

// Allocation-bearing convenience overload. Reuse an explicit plan when the
// same public shape is transformed repeatedly.
template <typename G, typename F, typename CoeffCtx>
void fft_forward_truncated_prefix(
        std::span<G> a,
        std::span<const F> roots,
        std::size_t xn,
        blst::TimingStats* timings = nullptr,
        std::string_view metric_prefix = "fft_trunc");

// Mutable only during construction. Callers that alter the public fields or
// schedules violate the executor precondition.
struct FftTruncatedForwardPlan
{
        struct StageWork
        {
                std::vector<std::size_t> left_j0;
                std::vector<std::size_t> left_nonzero;
                std::vector<std::size_t> twiddle_nonzero;
        };

        std::size_t n = 0;
        std::size_t xn = 0;
        std::size_t log2_n = 0;
        std::size_t butterflies_total = 0;
        std::size_t butterflies_executed = 0;
        std::size_t butterflies_skipped = 0;
        std::vector<StageWork> stage_work;
};

// Shape plan for an inverse transform that publishes only one output range.
// Mutable only during construction; concurrent execution requires read-only
// sharing of the completed plan.
struct FftInverseTruncatedRangePlan
{
        std::size_t n = 0;
        std::size_t start = 0;
        std::size_t length = 0;
        std::size_t stages = 0;
        std::size_t butterflies_total = 0;
        std::size_t butterflies_executed = 0;
        std::size_t butterflies_skipped = 0;
        std::vector<std::vector<std::size_t>> stage_left;
};

[[nodiscard]] constexpr std::size_t fft_log2_size(std::size_t n) noexcept
{
        std::size_t log2_n = 0;
        for (std::size_t t = n; t > 1; t >>= 1) {
                ++log2_n;
        }
        return log2_n;
}

inline FftTruncatedForwardPlan make_fft_truncated_forward_plan(std::size_t n, std::size_t xn)
{
        if (!blst::fft_detail::valid_transform_size(n) || xn == 0 || xn > n) {
                throw std::invalid_argument("invalid truncated forward FFT shape");
        }

        FftTruncatedForwardPlan plan;
        plan.n = n;
        plan.xn = xn;
        plan.log2_n = fft_log2_size(n);
        plan.butterflies_total = (n / 2) * plan.log2_n;

        if (xn == n) {
                plan.butterflies_executed = plan.butterflies_total;
                return plan;
        }

        plan.stage_work.resize(plan.log2_n);
        std::vector<unsigned char> active(n, 0);
        std::vector<unsigned char> next_active(n, 0);
        for (std::size_t i = 0; i < xn; ++i) {
                active[i] = 1;
        }
        bit_reverse<unsigned char>(std::span<unsigned char>(active));

        std::size_t stage_idx = 0;
        for (std::size_t len = 2; len <= n; len <<= 1, ++stage_idx) {
                const std::size_t half = len >> 1;
                const std::size_t stride = n / len;
                auto& stage = plan.stage_work[stage_idx];
                stage.left_j0.reserve(n / len);
                stage.left_nonzero.reserve((n / 2) - (n / len));
                stage.twiddle_nonzero.reserve((n / 2) - (n / len));
                std::fill(next_active.begin(), next_active.end(), static_cast<unsigned char>(0));

                for (std::size_t i = 0; i < n; i += len) {
                        for (std::size_t j = 0; j < half; ++j) {
                                const std::size_t left = i + j;
                                const std::size_t right = left + half;
                                if (!(active[left] || active[right])) {
                                        continue;
                                }
                                if (j == 0) {
                                        stage.left_j0.push_back(left);
                                }
                                else {
                                        stage.left_nonzero.push_back(left);
                                        stage.twiddle_nonzero.push_back(j * stride);
                                }
                                next_active[left] = 1;
                                next_active[right] = 1;
                        }
                }
                plan.butterflies_executed += stage.left_j0.size() + stage.left_nonzero.size();
                active.swap(next_active);
        }

        plan.butterflies_skipped = plan.butterflies_total - plan.butterflies_executed;
        return plan;
}

inline FftInverseTruncatedRangePlan make_fft_inverse_truncated_range_plan(
        std::size_t n,
        std::size_t start,
        std::size_t length)
{
        if (!blst::fft_detail::valid_transform_size(n) || length == 0 ||
            start >= n || length > n - start) {
                throw std::invalid_argument("invalid truncated inverse FFT range");
        }

        FftInverseTruncatedRangePlan plan;
        plan.n = n;
        plan.start = start;
        plan.length = length;
        plan.stages = fft_log2_size(n);
        plan.butterflies_total = (n / 2) * plan.stages;
        plan.stage_left.resize(plan.stages);

        std::vector<unsigned char> out_need(n, 0);
        std::vector<unsigned char> in_need(n, 0);
        for (std::size_t i = 0; i < length; ++i) {
                out_need[start + i] = 1;
        }

        for (std::size_t k = plan.stages; k > 0; --k) {
                const std::size_t len = std::size_t(1) << k;
                const std::size_t half = len >> 1;
                auto& lefts = plan.stage_left[k - 1];
                lefts.reserve(n / 2);
                std::fill(in_need.begin(), in_need.end(), static_cast<unsigned char>(0));

                for (std::size_t i = 0; i < n; i += len) {
                        for (std::size_t j = 0; j < half; ++j) {
                                const std::size_t left = i + j;
                                const std::size_t right = left + half;
                                if (!(out_need[left] || out_need[right])) {
                                        continue;
                                }
                                lefts.push_back(left);
                                in_need[left] = 1;
                                in_need[right] = 1;
                        }
                }

                plan.butterflies_executed += lefts.size();
                out_need.swap(in_need);
        }

        plan.butterflies_skipped = plan.butterflies_total - plan.butterflies_executed;
        return plan;
}

// Computes the cyclic length-n convolution by transforming A and, unless
// `skipBFft` is true, B in place, multiplying pointwise into C, and applying
// the inverse transform to C. All three buffers have exactly the same supported
// power-of-two length and are disjoint. Roots contain n forward powers and
// inv_roots contain n reciprocal powers. `xnA`/`xnB` are nonidentity prefix
// lengths; zero selects the full length, and every omitted tail is identity.
// With `skipBFft`, B must already be the complete transform of that padded
// prefix. On success C contains all n cyclic coefficients; A and B remain in
// the evaluation domain.
template <typename G1, typename G2, typename GT, typename F, typename CoeffCtx>
void poly_mult_fft(
	std::span<G1> A,
	std::span<G2> B,
	bool skipBFft,
	std::span<GT> C,
	std::span<const F> roots,
	std::span<const F> inv_roots,
	blst::TimingStats* timings = nullptr,
	std::size_t xnA = 0,
	std::size_t xnB = 0)
{
	std::size_t n = A.size();
	BLST_ASSERT(blst::fft_detail::valid_transform_size(n));
	BLST_ASSERT(B.size() == n);
	BLST_ASSERT(C.size() == n);
	BLST_ASSERT(roots.size() >= n);
	BLST_ASSERT(inv_roots.size() >= n);
	if (xnA == 0) {
		xnA = n;
	}
	if (xnB == 0) {
		xnB = n;
	}
	BLST_ASSERT(xnA > 0 && xnA <= n);
	BLST_ASSERT(xnB > 0 && xnB <= n);

	// Forward FFT on both polynomials (in-place)
	{
		blst::ScopedTimer timer(timings, "poly_mult.fft_a_ns");
		fft_forward_truncated_prefix<G1, F, CoeffCtx>(A, roots, xnA, timings, "poly_mult.fft_a");
	}

	if (skipBFft == false) {
		blst::ScopedTimer timer(timings, "poly_mult.fft_b_ns");
		fft_forward_truncated_prefix<G2, F, CoeffCtx>(B, roots, xnB, timings, "poly_mult.fft_b");
	}

	// Pointwise multiplication in evaluation domain
	{
		blst::ScopedTimer timer(timings, "poly_mult.pointwise_ns");
		for (std::size_t i = 0; i < n; ++i) {
			C[i] = CoeffCtx::mult(A[i], B[i]);
		}
	}

	// Inverse FFT to get product coefficients
	{
		blst::ScopedTimer timer(timings, "poly_mult.ifft_ns");
		fft<GT, F, CoeffCtx>(C, inv_roots, true);
	}
}

template <typename G, typename F, typename CoeffCtx>
void fft_forward_truncated_prefix(
        std::span<G> a,
        std::span<const F> roots,
        std::size_t xn,
        const FftTruncatedForwardPlan& plan,
        blst::TimingStats* timings,
        std::string_view metric_prefix) {
        const std::size_t n = a.size();
        BLST_ASSERT(blst::fft_detail::valid_transform_size(n));
        BLST_ASSERT(roots.size() >= n);
        BLST_ASSERT(xn > 0);
        BLST_ASSERT(xn <= n);
        BLST_ASSERT(plan.n == n);
        BLST_ASSERT(plan.xn == xn);
        BLST_ASSERT(plan.log2_n == fft_log2_size(n));
        BLST_ASSERT(xn == n || plan.stage_work.size() == plan.log2_n);

        if (xn == n) {
                fft<G, F, CoeffCtx>(a, roots, false);
                if (timings != nullptr) {
                        timings->add(metric_prefix, ".xn", xn);
                        timings->add(metric_prefix, ".butterflies_total", plan.butterflies_total);
                        timings->add(metric_prefix, ".butterflies_executed", plan.butterflies_executed);
                        timings->add(metric_prefix, ".butterflies_skipped", plan.butterflies_skipped);
                }
                return;
        }

        bit_reverse(a);

        for (std::size_t stage = 0, len = 2; len <= n; len <<= 1, ++stage) {
                const std::size_t half = len >> 1;
                const auto& stage_work = plan.stage_work[stage];

                for (std::size_t idx = 0; idx < stage_work.left_j0.size(); ++idx) {
                        const std::size_t left = stage_work.left_j0[idx];
                        const std::size_t right = left + half;
                        G u = a[left];
                        G t = a[right];
                        a[left] = CoeffCtx::plus(u, t);
                        a[right] = CoeffCtx::minus(u, t);
                }

                for (std::size_t idx = 0; idx < stage_work.left_nonzero.size(); ++idx) {
                        const std::size_t left = stage_work.left_nonzero[idx];
                        const std::size_t right = left + half;
                        const F& w = roots[stage_work.twiddle_nonzero[idx]];
                        G u = a[left];
                        G t = CoeffCtx::mult(a[right], w);
                        a[left] = CoeffCtx::plus(u, t);
                        a[right] = CoeffCtx::minus(u, t);
                }
        }

        if (timings != nullptr) {
                timings->add(metric_prefix, ".xn", xn);
                timings->add(metric_prefix, ".butterflies_total", plan.butterflies_total);
                timings->add(metric_prefix, ".butterflies_executed", plan.butterflies_executed);
                timings->add(metric_prefix, ".butterflies_skipped", plan.butterflies_skipped);
        }
}

template <typename G, typename F, typename CoeffCtx>
void fft_forward_truncated_prefix(
        std::span<G> a,
        std::span<const F> roots,
        std::size_t xn,
        blst::TimingStats* timings,
        std::string_view metric_prefix)
{
        const FftTruncatedForwardPlan plan = make_fft_truncated_forward_plan(a.size(), xn);
        fft_forward_truncated_prefix<G, F, CoeffCtx>(a, roots, xn, plan, timings, metric_prefix);
}

// Computes and normalizes only `[start, start + length)` of the inverse
// transform. Every frequency input is consumed. The requested range is valid
// on return; all other elements are unspecified scratch. The plan must match
// the exact public shape and may be shared only as immutable state.
template <typename G, typename F, typename CoeffCtx>
void fft_inverse_truncated_range(
        std::span<G> x,
        std::span<const F> inv_roots,
        std::size_t start,
        std::size_t length,
        const FftInverseTruncatedRangePlan& plan,
        blst::TimingStats* timings = nullptr,
        std::string_view metric_prefix = "ifft_trunc_range")
{
        const std::size_t n = x.size();
        BLST_ASSERT(blst::fft_detail::valid_transform_size(n));
        BLST_ASSERT(inv_roots.size() >= n);
        BLST_ASSERT(length > 0);
        BLST_ASSERT(start < n);
        BLST_ASSERT(length <= n - start);
        BLST_ASSERT(plan.n == n);
        BLST_ASSERT(plan.start == start);
        BLST_ASSERT(plan.length == length);
        BLST_ASSERT(plan.stages == fft_log2_size(n));
        BLST_ASSERT(plan.stage_left.size() == plan.stages);

        bit_reverse(x);

        for (std::size_t stage = 0, len = 2; len <= n; len <<= 1, ++stage) {
                const std::size_t half = len >> 1;
                const std::size_t stride = n / len;
                const auto& lefts = plan.stage_left[stage];

                for (std::size_t idx = 0; idx < lefts.size(); ++idx) {
                        const std::size_t left = lefts[idx];
                        const std::size_t right = left + half;
                        const std::size_t j = left & (half - 1);

                        if (j == 0) {
                                G u = x[left];
                                G t = x[right];
                                x[left] = CoeffCtx::plus(u, t);
                                x[right] = CoeffCtx::minus(u, t);
                        }
                        else {
                                const F& w = inv_roots[j * stride];
                                G u = x[left];
                                G t = CoeffCtx::mult(x[right], w);
                                x[left] = CoeffCtx::plus(u, t);
                                x[right] = CoeffCtx::minus(u, t);
                        }
                }
        }

        F inv_n = CoeffCtx::inverse(CoeffCtx::template from_u64<F>(n));
        for (std::size_t i = 0; i < length; ++i) {
                const std::size_t idx = start + i;
                x[idx] = CoeffCtx::mult(x[idx], inv_n);
        }

        if (timings != nullptr) {
                timings->add(metric_prefix, ".start", start);
                timings->add(metric_prefix, ".length", length);
                timings->add(metric_prefix, ".butterflies_total", plan.butterflies_total);
                timings->add(metric_prefix, ".butterflies_executed", plan.butterflies_executed);
                timings->add(metric_prefix, ".butterflies_skipped", plan.butterflies_skipped);
        }
}

// Allocation-bearing convenience overload. Reuse an explicit plan for
// repeated transforms of the same public output range.
template <typename G, typename F, typename CoeffCtx>
void fft_inverse_truncated_range(
        std::span<G> x,
        std::span<const F> inv_roots,
	std::size_t start,
	std::size_t length,
	blst::TimingStats* timings = nullptr,
	std::string_view metric_prefix = "ifft_trunc_range")
{
	const std::size_t n = x.size();
	BLST_ASSERT(blst::fft_detail::valid_transform_size(n));
	BLST_ASSERT(inv_roots.size() >= n);
	BLST_ASSERT(length > 0);
	BLST_ASSERT(start < n);
	BLST_ASSERT(length <= n - start);

	const FftInverseTruncatedRangePlan plan =
			make_fft_inverse_truncated_range_plan(n, start, length);
	fft_inverse_truncated_range<G, F, CoeffCtx>(
			x, inv_roots, start, length, plan, timings, metric_prefix);
}
