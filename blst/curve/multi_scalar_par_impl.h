#pragma once

#include "multi_scalar.h"

#ifdef BLST_CUDA
#include <cuda_runtime.h>

#ifndef __INTELLISENSE__
#include <thrust/iterator/iterator_facade.h>

#include <cub/cub.cuh>
//#include <cub/device/device_reduce.cuh>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>
#include <thrust/logical.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>

#endif
#include "blst/common/cuda_launch.h"
#endif

#include "blst/common/counting_iterator.h"
#include "blst/common/transform_iterator.h"
#include "blst/common/mdspan_compat.h"
#include <cstdint>
#include <vector>
#include <algorithm>

#include <inttypes.h>   // or <cinttypes>

namespace blst
{

	enum class Backend {
		CPU
#if defined(BLST_CUDA) && defined(BLST_CUDA_MSM)
		, CUDA
#endif
	};


	template<typename T>
	using d2Span = mdspan_compat::d2_span<T>;

	template<typename T>
	using d3Span = mdspan_compat::d3_span<T>;


	// cub will call operations on default-constructed Pointxyzz.
	// so we need to ensure that the default constructor is safe.
	// we also require alignment for Pointxyzz.
	template<typename Group>
	struct Defaulted : public Group {

		BLST_HD BLST_FORCEINLINE Defaulted()
		{
			memset(this, 0, sizeof(*this));
		}
		Defaulted(const Defaulted& other) = default;
		Defaulted& operator=(const Defaulted& other) = default;
	};

	BLST_HD BLST_INLINE auto msm_to_affine(GT& out, const GT& in)
	{
		out = in;
	}

	template<typename Field>
	BLST_HD BLST_INLINE auto msm_to_affine(AffinePoint<Field>& out, const Point<Field>& in)
	{
		Point<Field>::to_affine(out, in);
	}

	BLST_HD BLST_INLINE auto msm_from_affine(Defaulted<GT>& out, const GT& in)
	{
		static_cast<GT&>(out) = in;
	}
	template<typename Field>
	BLST_HD BLST_INLINE auto msm_from_affine(Defaulted<Pointxyzz<Field>>& out, const AffinePoint<Field>& in)
	{
		Pointxyzz<Field>::from_affine(out, in);
	}


	// ---- Packed val: (idx<<1)|sign --------------------------------------
	struct IdxSign {
		uint32_t v;
		// the index of the point in the input array that this value corresponds to
		BLST_HD BLST_FORCEINLINE uint32_t idx()  const { return v >> 1; }

		// the sign bit: 0 = positive, 1 = negative
		BLST_HD BLST_FORCEINLINE uint32_t sign() const { return v & 1; }
	};

#ifdef BLST_CUDA
	using Stream = cudaStream_t;
	using DeviceError = cudaError_t;
#else 
	using Stream = int;
	enum class DeviceError {
		cudaSuccess = 0,
		cudaErrorInvalidValue = 1
	};
#endif

	enum class MemcpyKind {
		HostToHost = 0, //Host->Host
		HostToDevice = 1, //Host->Device
		DeviceToHost = 2, //Device->Host
		DeviceToDevice = 3, //Device->Device
		Default = 4
	};

	template<Backend B>
	auto memcpyAsync(void* dst, const void* src, size_t count,
		MemcpyKind kind, Stream stream) -> DeviceError
	{
		if constexpr (B == Backend::CPU)
		{

			memcpy(dst, src, count);
			return DeviceError::cudaSuccess;
		}
		else
		{
#ifdef BLST_CUDA
			return ::cudaMemcpyAsync(dst, src, count, (cudaMemcpyKind)kind, stream);
#endif
		}
	}
	template<Backend B>
	auto memsetAsync(void* dst, int value, size_t count,
		Stream stream) -> DeviceError
	{
		if constexpr (B == Backend::CPU)
		{
			memset(dst, value, count);
			return DeviceError::cudaSuccess;
		}
		else
		{
#ifdef BLST_CUDA
			return ::cudaMemsetAsync(dst, value, count, stream);
#endif
		}
	}


	template<Backend B>
	auto streamSynchronize(Stream stream)
	{
		if constexpr (B == Backend::CPU)
		{
			return DeviceError::cudaSuccess;
		}
		else
		{
#ifdef BLST_CUDA
			return cudaStreamSynchronize(stream);
#endif // BLST_CUDA
		}
	}

	template<Backend B>
	auto getLastError()
	{
		if constexpr (B == Backend::CPU)
		{
			return DeviceError::cudaSuccess;
		}
		else
		{
#ifdef BLST_CUDA
			return cudaGetLastError();
#endif
		}
	}


	template<Backend B>
	auto getErrorString(DeviceError err)
	{
		if constexpr (B == Backend::CPU)
		{
			switch (err) {
			case DeviceError::cudaSuccess:
				return "cudaSuccess";
			case DeviceError::cudaErrorInvalidValue:
				return "cudaErrorInvalidValue";
			default:
				return "Unknown DeviceError";
			}
		}
		else
		{
#ifdef BLST_CUDA 
			return cudaGetErrorString(static_cast<cudaError_t>(err));
#endif
		}
	}

#define MSM_CHECK(call) do { \
	DeviceError err = call; \
	if (err != DeviceError::cudaSuccess) { \
		fprintf(stderr, "Error at %s:%d - %s\n", __FILE__, __LINE__, getErrorString<B>(err)); \
		BLST_ASSERT(0 && "MSM_CHECK failed"); \
	} }while(0)




	// Debug print function for digit extraction stage
	template<Backend B, typename Field>
	BLST_INLINE void print_extract_digit(
		const uint8_t* d_scalars,    // [M*N*nbytes] input scalars
		int nbytes,                  // bytes per scalar
		int nbits,                   // significant bits per scalar
		int wbits,                   // window size
		int N,                       // points per MSM
		int M,                       // number of MSMs
		int D,                       // number of windows
		d3Span<u8> d_digit,          // [D][M][N] extracted digits
		d3Span<IdxSign> d_index_signs, // [D][M][N] extracted index signs
		Stream stream)
	{
#ifndef NDEBUG
		printf("=== DIGIT EXTRACTION STAGE ===\n");
		printf("Parameters: M=%d, N=%d, D=%d, wbits=%d, nbits=%d, nbytes=%d\n", M, N, D, wbits, nbits, nbytes);

		// Allocate host memory for printing
		std::vector<uint8_t> h_scalars(M * N * nbytes);
		std::vector<u8> h_digits(D * M * N);
		std::vector<IdxSign> h_index_signs(D * M * N);

		// Copy data from device to host
		MSM_CHECK(memcpyAsync<B>(h_scalars.data(), d_scalars, M * N * nbytes,
			MemcpyKind::DeviceToHost, stream));
		MSM_CHECK(memcpyAsync<B>(h_digits.data(), d_digit.data_handle(), D * M * N * sizeof(u8),
			MemcpyKind::DeviceToHost, stream));
		MSM_CHECK(memcpyAsync<B>(h_index_signs.data(), d_index_signs.data_handle(), D * M * N * sizeof(IdxSign),
			MemcpyKind::DeviceToHost, stream));
		MSM_CHECK(streamSynchronize<B>(stream));

		// Print a sample of scalars
		printf("\nInput Scalars (first 16 per MSM):\n");
		for (int m = 0; m < M; m++) {
			printf("MSM %d:\n", m);
			for (int n = 0; n < N; n++) {
				printf("  Scalar[%d]: ", n);
				for (int b = 0; b < std::min(nbytes, 8); b++) {
					printf("%02x ", h_scalars[m * N * nbytes + n * nbytes + b]);
				}
				printf("\n");
			}
		}

		// Print extracted digits for first few elements
		printf("\nExtracted Digits:\n");
		for (int d = 0; d < D; d++) {
			for (int m = 0; m < M; m++) {
				printf("Window %d, MSM %d:\n", d, m);

				for (int n = 0; n < N; n++) {
					size_t global_idx = d * M * N + m * N + n;
					printf("  Point[%d]: digit=%u, idx=%u, sign=%u\n",
						n, h_digits[global_idx], h_index_signs[global_idx].idx(), h_index_signs[global_idx].sign());
				}
			}
		}

		//// Print digit distribution summary
		//std::vector<int> digit_counts(1 << wbits, 0);
		//for (size_t i = 0; i < h_digits.size(); i++) {
		//    if (h_digits[i] < digit_counts.size()) {
		//        digit_counts[h_digits[i]]++;
		//    }
		//}
		//printf("\nDigit distribution: ");
		//for (size_t i = 0; i < digit_counts.size() && i < 16; i++) {
		//    printf("d%zu:%d ", i, digit_counts[i]);
		//}
		printf("\n");
#endif
	}

	// Debug print function for sorting stage
	template<Backend B, typename Field>
	BLST_INLINE void print_sort_digit_indexsigns(
		int N,                       // points per MSM
		int M,                       // number of MSMs  
		int D,                       // number of windows
		int wbits,                   // window size
		d3Span<u8> d_digit,          // [D][M][N] original digits
		d3Span<IdxSign> d_index_signs, // [D][M][N] original index signs
		d3Span<u8> d_sorted_digit,   // [D][M][N] sorted digits
		d3Span<IdxSign> d_sorted_index_signs, // [D][M][N] sorted index signs
		Stream stream)
	{
#ifndef NDEBUG
		printf("=== SORTING STAGE ===\n");
		printf("Parameters: M=%d, N=%d, D=%d, wbits=%d\n", M, N, D, wbits);

		// Allocate host memory for printing
		std::vector<u8> h_sorted_digits(D * M * N);
		std::vector<IdxSign> h_sorted_index_signs(D * M * N);

		// Copy sorted data from device to host
		MSM_CHECK(memcpyAsync<B>(h_sorted_digits.data(), d_sorted_digit.data_handle(),
			D * M * N * sizeof(u8), MemcpyKind::DeviceToHost, stream));
		MSM_CHECK(memcpyAsync<B>(h_sorted_index_signs.data(), d_sorted_index_signs.data_handle(),
			D * M * N * sizeof(IdxSign), MemcpyKind::DeviceToHost, stream));
		MSM_CHECK(streamSynchronize<B>(stream));

		// Print sorted data for first segment
		printf("\nSorted data (first window, first MSM):\n");
		for (int d = 0; d < D; d++)
		{
			for (int m = 0; m < M; m++)
			{
				printf("Window %d, MSM %d:\n", d, m);
				size_t segment_offset = d * M * N + m * N;
				for (int i = 0; i < std::min(N, 10); i++) {
					size_t idx = segment_offset + i;
					printf("  [%d]: digit=%u, idx=%u, sign=%u\n",
						i, h_sorted_digits[idx], h_sorted_index_signs[idx].idx(), h_sorted_index_signs[idx].sign());
				}
			}
		}

		// Print sorting verification for each segment
		printf("\nSorting verification:\n");
		int sorted_segments = 0, total_segments = 0;
		for (int d = 0; d < D; d++) {
			for (int m = 0; m < M; m++) {
				total_segments++;
				size_t segment_offset = d * M * N + m * N;
				bool is_sorted = true;

				for (int i = 0; i < N - 1; i++) {
					if (h_sorted_digits[segment_offset + i] > h_sorted_digits[segment_offset + i + 1]) {
						is_sorted = false;
						break;
					}
				}

				if (is_sorted) sorted_segments++;
				printf("  Segment [d=%d][m=%d]: %s\n", d, m, is_sorted ? "SORTED" : "NOT SORTED");
			}
		}
		printf("Overall: %d/%d segments properly sorted\n", sorted_segments, total_segments);
#endif
	}

	// Debug print function for group reduction stage
	template<Backend B, typename Group, typename GroupAffine, typename GroupZZ>
	BLST_INLINE void print_reduce_groups(
		int N,                       // points per MSM
		int M,                       // number of MSMs  
		int D,                       // number of windows
		int W,                       // number of groups per MSM
		d3Span<u8> d_sorted_digit,   // [D][M][N] sorted digits
		d3Span<IdxSign> d_sorted_index_signs, // [D][M][N] sorted index signs
		d3Span<GroupZZ> d_groups,    // [D][M][W] reduced groups
		d2Span<const GroupAffine> d_points,   // [M*N] input points
		Stream stream)
	{
#ifndef NDEBUG
		printf("=== GROUP REDUCTION STAGE ===\n");
		printf("Parameters: M=%d, N=%d, D=%d, W=%d\n", M, N, D, W);

		// Allocate host memory for printing
		std::vector<u8> h_sorted_digits(D * M * N);
		std::vector<IdxSign> h_sorted_index_signs(D * M * N);
		std::vector<GroupZZ> h_groups(D * M * W);

		// Copy data from device to host
		MSM_CHECK(memcpyAsync<B>(h_sorted_digits.data(), d_sorted_digit.data_handle(),
			D * M * N * sizeof(u8), MemcpyKind::DeviceToHost, stream));
		MSM_CHECK(memcpyAsync<B>(h_sorted_index_signs.data(), d_sorted_index_signs.data_handle(),
			D * M * N * sizeof(IdxSign), MemcpyKind::DeviceToHost, stream));
		MSM_CHECK(memcpyAsync<B>(h_groups.data(), d_groups.data_handle(),
			D * M * W * sizeof(GroupZZ), MemcpyKind::DeviceToHost, stream));
		MSM_CHECK(streamSynchronize<B>(stream));

		// Print first few non-zero group coordinates
		printf("\ngroups [D][M][W]:\n");
		for (int d = 0; d < D; d++) {
			for (int m = 0; m < M; m++) {
				for (int g = 0; g < W; g++) {
					size_t group_idx = d * M * W + m * W + g;
					if (!vec_is_zero(&h_groups[group_idx], sizeof(h_groups[group_idx])))
					{
						Group point_affine;
						msm_to_jacobian(point_affine, h_groups[group_idx]);
						auto hash = hash_point<Group>(point_affine);
						printf("  Group[%d][%d][%d]: =%" PRIu64 "\n",
							d, m, g, hash);
					}
					else
					{
						printf("  Group[%d][%d][%d]: ZERO\n", d, m, g);
					}
				}
			}
		}
#endif
	}

	// Debug print function for bucket integration stage
	template<Backend B, typename Group, typename GroupZZ>
	BLST_INLINE void print_integrate_buckets_xyzz(
		int M,                       // number of MSMs
		int D,                       // number of windows  
		int W,                       // number of groups per MSM
		int wbits,                   // window size
		d3Span<GroupZZ> d_groups,    // [D][M][W] input groups
		d2Span<Group> d_window_sums,  // [D][M] output window sums
		Stream stream)
	{
#ifndef NDEBUG
		printf("=== BUCKET INTEGRATION STAGE ===\n");
		printf("Parameters: M=%d, D=%d, W=%d, wbits=%d\n", M, D, W, wbits);

		// Allocate host memory for printing
		std::vector<GroupZZ> h_groups(D * M * W);
		std::vector<Group> h_window_sums(D * M);

		// Copy data from device to host
		MSM_CHECK(memcpyAsync<B>(h_groups.data(), d_groups.data_handle(),
			D * M * W * sizeof(GroupZZ), MemcpyKind::DeviceToHost, stream));
		MSM_CHECK(memcpyAsync<B>(h_window_sums.data(), d_window_sums.data_handle(),
			D * M * sizeof(Group), MemcpyKind::DeviceToHost, stream));
		MSM_CHECK(streamSynchronize<B>(stream));

		// Print bucket integration input summary
		printf("\nInput buckets summary:\n");
		for (int d = 0; d < std::min(D, 2); d++) {
			for (int m = 0; m < std::min(M, 2); m++) {
				int non_zero = 0;
				size_t base = d * M * W + m * W;
				for (int w = 0; w < W; w++) {
					if (!vec_is_zero(&h_groups[base + w], sizeof(h_groups[base + w]))) {
						non_zero++;
					}
				}
				printf("  Window[%d][%d]: %d/%d non-zero buckets\n", d, m, non_zero, W);
			}
		}

		// Print window sums
		printf("\nWindow sums[D][M]:\n");
		for (int d = 0; d < D; d++) {
			for (int m = 0; m < M; m++) {
				size_t window_idx = d * M + m;
				bool is_identity = Group::is_identity(h_window_sums[window_idx]);
				if (!is_identity) {
					printf("  WindowSum[%d][%d]: %llu\n",
						d, m, (unsigned long long)hash_point(h_window_sums[window_idx]));
				}
				else {
					printf("  WindowSum[%d][%d]: identity\n", d, m);
				}
			}
		}

		// Statistics
		int identity_count = 0;
		int valid_count = 0;
		for (size_t i = 0; i < h_window_sums.size(); i++) {
			if (Group::is_identity(h_window_sums[i])) {
				identity_count++;
			}
			else {
				valid_count++;
			}
		}
		printf("\nWindow sum statistics: %d identity, %d valid\n", identity_count, valid_count);
#endif
	}

	// Debug print function for window folding stage
	template<Backend B, typename Group>
	BLST_INLINE void print_fold_windows(
		int M,                       // number of MSMs
		int D,                       // number of windows
		int wbits,                   // window size
		d2Span<Group> d_window_sums,  // [D][M] input window sums
		span<Group> d_outputs,     // [M] output results
		Stream stream)
	{
#ifndef NDEBUG
		printf("=== WINDOW FOLDING STAGE ===\n");
		printf("Parameters: M=%d, D=%d, wbits=%d\n", M, D, wbits);

		// Allocate host memory for printing
		std::vector<Group> h_window_sums(D * M);
		std::vector<Group> h_outputs(M);

		// Copy data from device to host
		MSM_CHECK(memcpyAsync<B>(h_window_sums.data(), d_window_sums.data_handle(),
			D * M * sizeof(Group), MemcpyKind::DeviceToHost, stream));
		MSM_CHECK(memcpyAsync<B>(h_outputs.data(), d_outputs.data(),
			M * sizeof(Group), MemcpyKind::DeviceToHost, stream));
		MSM_CHECK(streamSynchronize<B>(stream));

		// Print input window sums for folding
		printf("\nInput window sums for folding[M][D]:\n");
		for (int m = 0; m < M; m++) {
			printf("  MSM %d:\n", m);
			for (int d = 0; d < D; d++) {
				size_t idx = d * M + m;
				bool is_identity = Group::is_identity(h_window_sums[idx]);
				if (!is_identity) {
					printf("    Window[%d]: %llu\n",
						d, (unsigned long long)hash_point(h_window_sums[idx]));
				}
				else {
					printf("    Window[%d]: identity\n", d);
				}
			}
		}

		// Print final outputs
		printf("\nFinal MSM outputs[M]:\n");
		for (int m = 0; m < M; m++) {
			bool is_identity = Group::is_identity(h_outputs[m]);

			if (!is_identity) {
				printf("  MSM[%d]: %llu\n",
					m, (unsigned long long)hash_point(h_outputs[m]));
			}
			else {
				printf("  MSM[%d]: identity\n", m);
			}
		}
#endif
	}



























	// ========== Kernels ==========
	BLST_HD BLST_INLINE void extract_digits_flat(
		int window_idx,
		int i,
		const uint8_t* BLST_restrict scalars,
		int N,                // number of scalars
		int nbytes,           // bytes per scalar
		int nbits,            // significant bits per scalar
		WindowExtract* d_windows,
		int window_size,            // window size (e.g., 5). W=2^window_size
		int D,                // number of windows = ceil(nbits / window_size)
		u8* BLST_restrict digits,
		IdxSign* BLST_restrict index_signs)
	{
		BLST_ASSERT(window_idx * window_size <= nbits && "the start of the window can at most be the bit after then last");
		auto outIdx = window_idx * N + i;

		const uint8_t* s = scalars + i * nbytes;
		WindowExtract& we = d_windows[window_idx];
		auto encoding = we.extract(s);

		u8 digit = encoding & ((limb_t)1 << we.encode_size) - 1;
		u8 sign = (encoding >> we.encode_size) & 1;

		digits[outIdx] = digit;
		index_signs[outIdx] = IdxSign{ u32(sign) | (i << 1) };
	}

#ifdef BLST_CUDA
	// taka flat array of scalars, extract the w-bit Booth encoded digits
	// the output order is column-major by digit. First row has the N LSB digits, ...
	BLST_GLOBAL inline void extract_digits_flat_kernel(
			const uint8_t* BLST_restrict scalars,
			int N,                // number of scalars
			int nbytes,           // bytes per scalar
			int nbits,            // significant bits per scalar
			WindowExtract* d_windows,
			int wbits,            // window size (e.g., 5). W=2^wbits
			int D,                // number of windows = ceil(nbits / wbits)
			u8* BLST_restrict digits,
			IdxSign* BLST_restrict index_signs)   // [0,N]|{0,1}^N*D
	{
		int tid = blockDim.x * blockIdx.x + threadIdx.x;
		int total = N * D;
		if (tid >= total)
			return;

		int i = tid / D;     // scalar index
		int d = tid % D;     // window index (0 = LSB window)

		extract_digits_flat(d, i, scalars, N, nbytes, nbits, d_windows, wbits, D, digits, index_signs);
	}

#endif // BLST_CUDA

	// taka flat array of scalars, extract the w-bit Booth encoded digits
	// the output order is column-major by digit. First row has the N LSB digits, ...
	template<Backend B>
	DeviceError extract_digits_flat(
		d2Span<const uint8_t> scalars, // [M*N][nbytes]
		int num_scalars,                // number of scalars
		span<WindowExtract> d_windows, // [D]
		int nbytes,           // bytes per scalar
		int nbits,            // significant bits per scalar
		int wbits,            // window size (e.g., 5). W=2^wbits
		int D,                // number of windows = ceil(nbits / wbits)
		d2Span<u8> digit,    // [D][M*N]
		d2Span<IdxSign> idx_sign, // [D][M*N]
		Stream stream)
	{
		if (num_scalars <= 0 || D <= 0)
			return DeviceError::cudaSuccess;

		BLST_ASSERT(wbits > 0 && wbits <= 8);
		BLST_ASSERT(d_windows.size() == D);
		BLST_ASSERT(digit.extent(0) == static_cast<size_t>(D));
		BLST_ASSERT(idx_sign.extent(0) == static_cast<size_t>(D));
		BLST_ASSERT(digit.extent(1) == static_cast<size_t>(num_scalars));
		BLST_ASSERT(idx_sign.extent(1) == static_cast<size_t>(num_scalars));


		int total = num_scalars * D;
		int blockSize = 256;
		int numBlocks = (total + blockSize - 1) / blockSize;

		if constexpr (B == Backend::CPU)
		{
			for (int tid = 0; tid < total; tid++)
			{
				int i = tid / D;     // scalar index
				int d = tid % D;     // window index (0 = LSB window)
				extract_digits_flat(
					d, i,
					scalars.data_handle(), num_scalars,
					nbytes, nbits,
					d_windows.data(), wbits, D,
					digit.data_handle(),
					idx_sign.data_handle());
			}
			return DeviceError::cudaSuccess;
		}
		else
		{
#ifdef BLST_CUDA
			extract_digits_flat_kernel BLST_LAUNCH(numBlocks, blockSize, 0, stream) (
				scalars.data_handle(), num_scalars,
				nbytes, nbits,
				d_windows.data(), wbits, D,
				digit.data_handle(),
				idx_sign.data_handle());
#endif // BLST_CUDA
		}


		return getLastError<B>();
	}



	template< Backend B,
		typename KeyT, typename ValueT, typename BeginOffsetIteratorT, typename EndOffsetIteratorT>
	static inline DeviceError SortPairs(
		void* d_temp_storage,
		size_t& temp_storage_bytes,
		const KeyT* d_keys_in,
		KeyT* d_keys_out,
		const ValueT* d_values_in,
		ValueT* d_values_out,
		u64 num_items,
		u64 num_segments,
		BeginOffsetIteratorT d_begin_offsets,
		EndOffsetIteratorT d_end_offsets,
		int begin_bit = 0,
		int end_bit = sizeof(KeyT) * 8,
		Stream stream = 0)
	{
		if constexpr (B == Backend::CPU)
		{
			if (d_temp_storage == nullptr)
			{
				temp_storage_bytes = 1; // dummy value
				return DeviceError::cudaSuccess;
			}


			//while (d_begin_offsets != d_end_offsets)
			for (u64 segment = 0; segment != num_segments; ++segment)
			{
				auto begin = *d_begin_offsets++;
				auto end = *d_end_offsets++;

				std::vector<std::pair<KeyT, ValueT>> items;
				for (u64 i = begin; i < end; i++)
				{
					items.push_back({ d_keys_in[i],d_values_in[i] });
				}
				std::sort(items.begin(), items.end(),
					[=](const std::pair<KeyT, ValueT>& a, const std::pair<KeyT, ValueT>& b)
					{
						return a.first < b.first;
					});
				for (u64 i = begin; i < end; i++)
				{
					d_keys_out[i] = items[i - begin].first;
					d_values_out[i] = items[i - begin].second;
				}
			}
			return DeviceError::cudaSuccess;
		}
		else
		{
#ifdef BLST_CUDA
			return cub::DeviceSegmentedRadixSort::SortPairs(
				d_temp_storage,
				temp_storage_bytes,
				d_keys_in,
				d_keys_out,
				d_values_in,
				d_values_out,
				num_items,
				num_segments,
				d_begin_offsets,
				d_end_offsets,
				begin_bit,
				end_bit,
				stream);
#endif // BLST_CUDA
		}
	}


	// Sort each of D digit slices of N elements with a single segmented sort.
	// given D digit arrays of N elements each, sort each digit array
	// we will sort key value pairs, where the key is the digit value (u8)
	// The value is the index of the point in the input array that this digit
	// corresponds to. The LSB of the value is reserved for the sign bit.
	//
	template<Backend B>
	DeviceError sort_digits_segmented(
		void* d_workspace,
		size_t& workspace_bytes,
		const u8* d_keys_input,          // [D*N], column-major by digit
		u8* d_keys_output,          // [D*N], column-major by digit
		const IdxSign* d_vals_input,          // [D*N]
		IdxSign* d_vals_output,          // [D*N]
		int          N,
		int          D,
		int          wbits,
		Stream stream = 0)
	{
		if (N <= 0 || D <= 0) return DeviceError::cudaSuccess;
		if (wbits < 0) return DeviceError::cudaErrorInvalidValue;
		if (wbits > 8) return DeviceError::cudaErrorInvalidValue;

		size_t count = static_cast<size_t>(N) * static_cast<size_t>(D);
		auto d_begins = transform_iterator(
			counting_iterator<u32>(0),
			[=] BLST_HD(u32 i) { return i * N; }
		);
		auto d_ends = d_begins + 1;


		// 1) Query CUB temp storage
		size_t required = 0;
		SortPairs<B>(
			nullptr,                           // void *d_temp_storage,
			required,                           // size_t &temp_storage_bytes,
			d_keys_input,                       // const KeyT *d_keys_in,
			d_keys_output,                       // KeyT *d_keys_out,
			d_vals_input,                       // const ValueT *d_values_in,
			d_vals_output,                       // ValueT *d_values_out,
			count,                               // ::cuda::std::int64_t num_items,
			D,                                   // ::cuda::std::int64_t num_segments,
			d_begins,                           // BeginOffsetIteratorT d_begin_offsets,
			d_ends,                            // EndOffsetIteratorT d_end_offsets,
			0,                                   // int begin_bit = 0,
			wbits,                               // int end_bit = sizeof(KeyT) * 8,
			stream                               // Stream stream = 0,
		);

		if (d_workspace == nullptr) {
			workspace_bytes = required;
			return DeviceError::cudaSuccess;
		}
		if (workspace_bytes < required)
			return DeviceError::cudaErrorInvalidValue;

		// 2) Run segmented sort
		auto ec = SortPairs<B>(
			d_workspace,
			workspace_bytes,
			d_keys_input,
			d_keys_output,
			d_vals_input,
			d_vals_output,
			count,
			D,
			d_begins,
			d_ends,
			0,
			wbits,
			stream);
		if (ec != DeviceError::cudaSuccess)
			return ec;

		return DeviceError::cudaSuccess;
	}



	// ---- Transform functor: packed -> Pointxyzz by gathering -------------
	template<class GroupAffine, typename GroupZZ>
	struct PackedToPointXYZZ {
		const GroupAffine* BLST_restrict points; // device
		int num_points;

		BLST_HD BLST_FORCEINLINE
			GroupZZ operator()(IdxSign index_sign) const {

			BLST_ASSERT(index_sign.idx() < num_points);
			GroupAffine p = points[index_sign.idx()];
			GroupAffine::cneg(p, index_sign.sign());
			GroupZZ q;
			msm_from_affine(q, p);
			return q;
		}
	};

	// ---- Point monoid for reduces/scans ----------------------------------
	// add: DefaultedPointxyzz x DefaultedPointxyzz -> DefaultedPointxyzz
	template<class Group>
	struct XYZZAdd {
		BLST_HD BLST_FORCEINLINE
			Group operator()(
				auto&& a,
				auto&& b) const
		{
			Group r;
			Group::dadd(r, a, b);
			return r;
		}
	};

	// Helpers: align up sizes / pointers
	BLST_HD BLST_INLINE static size_t align_up_size(
		size_t x,
		size_t a)
	{
		size_t mask = a - 1;
		return (x + mask) & ~mask;
	}

	BLST_HD BLST_INLINE static uint8_t* align_up_ptr(
		uint8_t* p,
		size_t a)
	{
		uintptr_t u = reinterpret_cast<uintptr_t>(p);
		uintptr_t aligned = (u + (a - 1)) & ~(uintptr_t)(a - 1);
		return reinterpret_cast<uint8_t*>(aligned);
	}


	template<typename GroupZZ, typename ValT>
	BLST_HD BLST_INLINE void scatter_group_sums_multi_kernel(
		u64 tid,
		const uint32_t* __restrict d_group_keys,
		const ValT* __restrict d_group_vals,
		const size_t* __restrict d_num_runs_ptr,
		int W,
		int num_sets,
		GroupZZ* __restrict d_groups)
	{
		if (tid >= *d_num_runs_ptr)
			return;

		uint32_t gk = d_group_keys[tid];
		uint32_t wspan = static_cast<uint32_t>(W + 1);
		uint32_t m = gk / wspan;
		uint32_t k = gk % wspan;

		if (k == 0)
			return;

		size_t out_idx = static_cast<size_t>(m) * static_cast<size_t>(W) + static_cast<size_t>(k - 1);
		//if (out_idx >= num_sets * W)
		//{
		//	printf("ERROR: scatter_group_sums_multi_kernel invalid output index %zu (m=%u, k=%u, W=%d, M=%d)\n",
		//		out_idx, m, k, W, num_sets);
		//	return;
		//}
		d_groups[out_idx] = static_cast<GroupZZ>(d_group_vals[tid]);
	}


#ifdef BLST_CUDA


	template<typename GroupZZ, typename ValT>
	__global__ void scatter_group_sums_multi_kernel(
		const uint32_t* __restrict d_group_keys,
		const ValT* __restrict d_group_vals,
		const size_t* __restrict d_num_runs_ptr,
		int W,
		int num_sets,
		GroupZZ* __restrict d_groups)
	{
		size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
		scatter_group_sums_multi_kernel<GroupZZ, ValT>(
			tid,
			d_group_keys,
			d_group_vals,
			d_num_runs_ptr,
			W,
			num_sets,
			d_groups);
	}
#endif // BLST_CUDA

	// copy each group to its output location in d_groups
	// d_group_keys: [num_sets * (W+1)] in [num_sets] * {0,1,...,W} , in packed form.
	// d_group_vals: [num_sets * (W+1)]
	// d_num_runs: [1], the number of actual groups, excluding empty groups

	template<Backend B,
		typename GroupZZ, typename ValT>
	DeviceError scatter_group_sums_multi(
		span<const uint32_t> d_group_keys,   // [num_sets * (W+1)], d_num_runs of them are valid
		span<const ValT> d_group_vals,       // [num_sets * (W+1)], d_num_runs of them are valid
		const size_t* __restrict d_num_runs, // [1], the number of actual groups, excluding empty groups
		int W,                               // the number of groups per set (max)
		int num_sets,                        // the number of sets (D*M)
		d2Span<GroupZZ> d_groups, // [num_sets][W], zero groups skipped.
		Stream stream)
	{
		if constexpr (B == Backend::CPU)
		{
			for (size_t tid = 0; tid < *d_num_runs; tid++)
			{
				scatter_group_sums_multi_kernel<GroupZZ, ValT>(
					tid,
					d_group_keys.data(),
					d_group_vals.data(),
					d_num_runs,
					W,
					num_sets,
					d_groups.data_handle());
			}
			return DeviceError::cudaSuccess;

		}
		else
		{
#if defined(BLST_CUDA) && defined(BLST_CUDA_MSM)
			static_assert(B == Backend::CUDA);
			int blk = 256;
			size_t scatter_threads = W;
			int grd = static_cast<int>((scatter_threads + blk - 1) / blk);

			scatter_group_sums_multi_kernel<GroupZZ, ValT>
				BLST_LAUNCH(grd, blk, 0, stream) (
					d_group_keys.data(),
					d_group_vals.data(),
					d_num_runs,
					W,
					num_sets,
					d_groups.data_handle());

			return getLastError<B>();
#endif // BLST_CUDA
		}

	}

	template<Backend B,
		typename KeysInputIteratorT, typename UniqueOutputIteratorT, typename ValuesInputIteratorT, typename AggregatesOutputIteratorT, typename NumRunsOutputIteratorT, typename ReductionOpT, typename NumItemsT>
	static inline DeviceError ReduceByKey(
		void* d_temp_storage,
		size_t& temp_storage_bytes,
		KeysInputIteratorT d_keys_in,
		UniqueOutputIteratorT d_unique_out,
		ValuesInputIteratorT d_values_in,
		AggregatesOutputIteratorT d_aggregates_out,
		NumRunsOutputIteratorT d_num_runs_out,
		ReductionOpT reduction_op,
		NumItemsT num_items,
		Stream stream = 0)
	{
		if constexpr (B == Backend::CPU)
		{

			if (!d_temp_storage)
			{
				temp_storage_bytes = 1; // dummy value
				return DeviceError::cudaSuccess;
			}

			// Simple CPU implementation of ReduceByKey
			using KeyT = typename std::iterator_traits<KeysInputIteratorT>::value_type;
			using ValT = typename std::iterator_traits<ValuesInputIteratorT>::value_type;
			if (num_items == 0)
			{
				if (d_num_runs_out != nullptr)
				{
					*d_num_runs_out = 0;
				}
				return DeviceError::cudaSuccess;
			}
			KeyT current_key = *(d_keys_in);
			ValT current_aggregate = *(d_values_in);
			size_t run_count = 0;
			for (NumItemsT i = 1; i < num_items; i++)
			{
				KeyT next_key = *(d_keys_in + i);
				ValT next_value = *(d_values_in + i);
				if (next_key == current_key)
				{
					current_aggregate = reduction_op(current_aggregate, next_value);
				}
				else
				{
					*(d_unique_out + run_count) = current_key;
					*(d_aggregates_out + run_count) = current_aggregate;
					run_count++;
					current_key = next_key;
					current_aggregate = next_value;
				}
			}
			// Write the last run
			*(d_unique_out + run_count) = current_key;
			*(d_aggregates_out + run_count) = current_aggregate;
			run_count++;
			if (d_num_runs_out != nullptr)
			{
				*d_num_runs_out = run_count;
			}
			return DeviceError::cudaSuccess;
		}
		else
		{
#if defined(BLST_CUDA) && defined(BLST_CUDA_MSM)
			static_assert(B == Backend::CUDA);
			return cub::DeviceReduce::ReduceByKey(
				d_temp_storage,
				temp_storage_bytes,
				d_keys_in,
				d_unique_out,
				d_values_in,
				d_aggregates_out,
				d_num_runs_out,
				reduction_op,
				num_items,
				stream);
#endif
		}

	}

	template< Backend B,
		typename InputIterator, typename Predicate>
	bool all_of(
		Stream stream,
		InputIterator first,
		InputIterator last,
		Predicate pred)
	{
		if constexpr (B == Backend::CPU)
		{
			for (auto it = first; it != last; it++)
			{
				if (!pred(*it))
					return false;
			}
			return true;
		}
		else
		{
#if defined(BLST_CUDA) && defined(BLST_CUDA_MSM)

			static_assert(B == Backend::CUDA);
			return thrust::all_of(
				thrust::cuda::par.on(stream),
				first,
				last,
				pred);
#endif
		}
	}


	// ===== Per-group reduce: keys/vals sorted by key =======================
	// 
	// Takes an input M vectors of size N, each with a key in [0..W].
	// Groups are defined by keys. Key 0 is skipped.
	// For each M, output W groups, each being the sum of all points
	// within the respective group.
	// 
	// @d_workspace: temp buffer, or nullptr to query size c
	// @workspace_bytes: in: size of d_workspace; out: required size if d_workspace is nullptr or too small
	// @d_keys: size M * N matrix iter with values in [0..W], in grouped order
	// @d_vals: size M * N matrix iter of points to be reduced.
	// @d_groups: size M * W vec of Group, the output groups, skipping group 0.
	// @M: number of independent vectors (digit vectors)
	// @N: number of keys/vals/points
	// @W: number of groups = 2^(wbits-1). keys are in [0..W], W+1 possible values. 0 is skipped.
	// @stream: cuda stream
	// 
	// This algorithm works as follows:
	//   1) load the points: we have N vals which are packed (idx<<1)|sign,
	//      for each val, we load points[idx] and apply sign to get DefaultedPointxyzz.
	//   2) for each group (key), we reduce the points in that group to get 
	//      the group sum. This will be done via  thrust::reduce_by_key.
	//   3) We scatter these to the output groups array. we will ignore group 0 
	//      and fill any missing groups with the point at infinity.
	//
//	template<class GroupZZ, typename SrcIter>
//	DeviceError reduce_groups_global(
//		void* d_workspace,
//		size_t& workspace_bytes,
//		d2Span<const u8> d_keys, // [D*M][N] keys in [0..W]
//		SrcIter&& src_iter,      // flat [D*M][N] input iterater, *iter = GroupZZ
//		d2Span<GroupZZ> d_groups,// [D*M][W] output groups
//		int num_sets, // D*M
//		int N,
//		int W,
//		Stream stream)
//	{
//
//
//		using SrcItDecayed = typename std::decay<SrcIter>::type;
//		using ValT = ::blst::Defaulted<GroupZZ>;
//		using deref_type = typename std::iterator_traits<SrcItDecayed>::value_type;
//		static_assert(std::is_same<deref_type, ValT>::value,
//			"reduce_groups: src_iter must deref to Defaulted<GroupZZ>");
//
//		size_t total = num_sets * N;
//
//		// the number of groups, including group 0
//		size_t total_groups = num_sets * (W + 1);
//
//
//		auto idx_first = counting_iterator<size_t>(0);
//		auto global_key_iter = transform_iterator(
//			idx_first,
//			[=] BLST_HD(size_t idx) -> uint32_t
//		{
//			u8 k = d_keys.data_handle()[idx];
//			size_t m = idx / static_cast<size_t>(N);
//			size_t wspan = static_cast<size_t>(W + 1);
//			return static_cast<uint32_t>(m * wspan + static_cast<size_t>(k));
//		});
//
//		size_t temp_storage_bytes = 0;
//		{
//			uint32_t* dummy_keys = nullptr;
//			ValT* dummy_vals = nullptr;
//			size_t* dummy_num_runs = nullptr;
//
//			DeviceError err = ReduceByKey(
//				nullptr,                   // void *d_temp_storage,
//				temp_storage_bytes,        // size_t &temp_storage_bytes,
//				global_key_iter,           // KeysInputIteratorT d_keys_in,
//				dummy_keys,                // UniqueOutputIteratorT d_unique_out,
//				src_iter,                  // ValuesInputIteratorT d_values_in,
//				dummy_vals,                // AggregatesOutputIteratorT d_aggregates_out,
//				dummy_num_runs,            // NumRunsOutputIteratorT d_num_runs_out,
//				XYZZAdd<ValT>{},          // ReductionOpT reduction_op,
//				static_cast<int>(total),   // NumItemsT num_items,
//				stream);                   // Stream stream = 0,
//			if (err != DeviceError::cudaSuccess)
//				return err;
//
//			size_t required = 0;
//			required += temp_storage_bytes;
//			required = align_up_size(required, alignof(ValT));
//			required += sizeof(ValT) * total_groups;
//			required = align_up_size(required, alignof(uint32_t));
//			required += sizeof(uint32_t) * total_groups;
//			required = align_up_size(required, alignof(size_t));
//			required += sizeof(size_t);
//
//			if (d_workspace == nullptr)
//			{
//				workspace_bytes = required;
//				return DeviceError::cudaSuccess;
//			}
//
//			if (workspace_bytes < required)
//			{
//				workspace_bytes = required;
//				return DeviceError::cudaErrorInvalidValue;
//			}
//		}
//
//
//
//#ifndef NDEBUG
//		{
//			bool keys_ok = all_of(stream,
//				d_keys.data_handle(),
//				d_keys.data_handle() + total,
//				[=] BLST_HD(u8 k)
//			{
//				return k <= static_cast<u8>(W);
//			});
//			if (!keys_ok)
//				return DeviceError::cudaErrorInvalidValue;
//		}
//#endif
//		uint8_t* ws = static_cast<uint8_t*>(d_workspace);
//
//		void* d_temp_storage = ws;
//		ws += temp_storage_bytes;
//
//		ws = align_up_ptr(ws, alignof(ValT));
//
//		// flat buffer for group values, some might be missing
//		span<ValT> d_group_vals(reinterpret_cast<ValT*>(ws), num_sets * (W + 1)); // [D*M * (W+1)]
//		ws += sizeof(ValT) * total_groups;
//
//		ws = align_up_ptr(ws, alignof(uint32_t));
//		// flat buffer for group keys, some might be missing
//		span<uint32_t> d_group_keys(reinterpret_cast<uint32_t*>(ws), num_sets * (W + 1)); // [D*M * (W+1)]
//		ws += sizeof(uint32_t) * total_groups;
//
//		// buffer to contain the actual number of groups found
//		ws = align_up_ptr(ws, alignof(size_t));
//		size_t* d_num_runs = reinterpret_cast<size_t*>(ws);
//		ws += sizeof(size_t);
//
//		BLST_ASSERT(num_sets > 0);
//		BLST_ASSERT(d_keys.extent(0) == num_sets);
//		BLST_ASSERT(d_keys.extent(1) == N);
//		BLST_ASSERT(d_groups.extent(0) == num_sets);
//		BLST_ASSERT(d_groups.extent(1) == W);
//
//
//
//		DeviceError err = ReduceByKey(
//			d_temp_storage,
//			temp_storage_bytes,
//			global_key_iter,
//			d_group_keys.data(),
//			src_iter,
//			d_group_vals.data(),
//			d_num_runs,
//			XYZZAdd<ValT>{},
//			static_cast<int>(total),
//			stream);
//		if (err != DeviceError::cudaSuccess)
//			return err;
//
//		//for (u64 i = 0; i < *d_num_runs; i++)
//		//{
//		//	if (d_group_keys[i] >= total_groups)
//		//	{
//		//		printf("ERROR: reduce_groups invalid group key %u (total_groups=%zu)\n",
//		//			d_group_keys[i], total_groups);
//		//		return DeviceError::cudaErrorInvalidValue;
//		//	}
//		//	Group p;
//		//	GroupZZ::to_Jacobian(p, d_group_vals[i]);
//		//	printf("  group key %u -> %llu\n",
//		//		d_group_keys[i], hash_point(p));
//		//}
//
//		// since some groups may be missing, we need to scatter the results
//		// and zero out the missing groups.
//		size_t out_elems = static_cast<size_t>(num_sets) * static_cast<size_t>(W);
//		auto ec = memsetAsync(d_groups.data_handle(), 0, sizeof(GroupZZ) * out_elems, stream);
//		if (ec != DeviceError::cudaSuccess)
//			return ec;
//
//
//		// scatter the group sums to their output locations, skipping group 0
//		// and skipping missing groups.
//		ec = scatter_group_sums_multi<GroupZZ, ValT>(
//			d_group_keys,
//			d_group_vals,
//			d_num_runs,
//			W,
//			num_sets,
//			d_groups,
//			stream);
//
//		return ec;
//	}

	template<Backend B,
		typename InputIterator1, typename InputIterator2, typename OutputIterator1, typename OutputIterator2, typename BinaryPredicate, typename BinaryFunction>
	std::pair<OutputIterator1, OutputIterator2> reduce_by_key(
		Stream stream,
		InputIterator1 keys_first,
		InputIterator1 keys_last,
		InputIterator2 values_first,
		OutputIterator1 keys_output,
		OutputIterator2 values_output,
		BinaryPredicate binary_pred,
		BinaryFunction binary_op
	)
	{
		if constexpr (B == Backend::CPU)
		{

			// cpu implementation of reduce_by_key
			u8* dummy = (u8*)1;
			size_t temp_storage_size = 0;
			u64 d_num_runs_out;
			ReduceByKey<B>(dummy, temp_storage_size,
				keys_first,
				keys_output,
				values_first,
				values_output,
				&d_num_runs_out,
				binary_op,
				keys_last - keys_first,
				stream);
			return { keys_output + d_num_runs_out, values_output + d_num_runs_out };

		}
		else
		{
#if defined(BLST_CUDA) && defined(BLST_CUDA_MSM)

			static_assert(B == Backend::CUDA);

			auto policy = thrust::cuda::par.on(stream);
			return thrust::reduce_by_key(
				policy,
				keys_first,
				keys_last,
				values_first,
				keys_output,
				values_output,
				binary_pred,
				binary_op);
#endif
		}
	}

	struct Equal {
		BLST_HD BLST_FORCEINLINE
			bool operator()(u8 a, u8 b) const {
			return a == b;
		}
	};

	template<Backend B,
		class GroupZZ, typename SrcIter>
	DeviceError reduce_groups(
		void* d_workspace,
		size_t& workspace_bytes,
		d2Span<const u8>           d_keys,          // [num_sets][N], keys in [0..W]
		SrcIter&& src_iter,        // flat [num_sets * N], *it = Defaulted<GroupZZ>
		d2Span<GroupZZ>            d_groups,        // [num_sets][W]
		int                        num_sets,
		int                        N,
		int                        W,
		Stream                     stream)
	{
		using ValT = ::blst::Defaulted<GroupZZ>;

		// Carve workspace for temporary buffers
		uint8_t* ws_ptr = static_cast<uint8_t*>(d_workspace);

		// Align and allocate unique keys buffer
		ws_ptr = align_up_ptr(ws_ptr, alignof(u32));
		u32* uniq_keys = reinterpret_cast<u32*>(ws_ptr);
		ws_ptr += sizeof(u32) * (W + 1);

		// Align and allocate aggregated values buffer
		ws_ptr = align_up_ptr(ws_ptr, alignof(ValT));
		ValT* agg_vals = reinterpret_cast<ValT*>(ws_ptr);
		ws_ptr += sizeof(ValT) * (W + 1);

		size_t* d_runs;
		// Align and allocate runs buffer
		ws_ptr = align_up_ptr(ws_ptr, alignof(size_t));
		d_runs = reinterpret_cast<size_t*>(ws_ptr);
		ws_ptr += sizeof(size_t);


		// Ensure we haven't exceeded workspace bounds
		auto ws_needed = static_cast<size_t>(ws_ptr - static_cast<uint8_t*>(d_workspace));

		if (d_workspace == nullptr)
		{
			workspace_bytes = std::max(workspace_bytes, ws_needed);
			return DeviceError::cudaSuccess;
		}
		else
		{
			if (ws_needed > workspace_bytes)
			{
				workspace_bytes = ws_needed;
				return DeviceError::cudaErrorInvalidValue;
			}
		}

		BLST_ASSERT(d_keys.extent(0) == (size_t)num_sets);
		BLST_ASSERT(d_keys.extent(1) == (size_t)N);
		BLST_ASSERT(d_groups.extent(0) == (size_t)num_sets);
		BLST_ASSERT(d_groups.extent(1) == (size_t)W);

		size_t total_out = (size_t)num_sets * (size_t)W;

		MSM_CHECK(memsetAsync<B>(
			d_groups.data_handle(),
			0,
			sizeof(GroupZZ) * total_out,
			stream));


		for (int s = 0; s < num_sets; s++)
		{
			const u8* keys_begin = d_keys.data_handle() + (size_t)s * (size_t)N;
			const u8* keys_end = keys_begin + N;

			auto vals_begin = src_iter + (size_t)s * (size_t)N;

			auto end_pair = reduce_by_key<B>(
				stream,				 		// Stream stream,
				keys_begin,			 		// InputIterator1 keys_first,
				keys_end,			 		// InputIterator1 keys_last,
				vals_begin,			 		// InputIterator2 values_first,
				uniq_keys,			 		// OutputIterator1 keys_output,
				agg_vals,			 		// OutputIterator2 values_output,
				Equal{},		 				// BinaryPredicate binary_pred,
				XYZZAdd<ValT>{});	 		// BinaryFunction binary_op

			size_t runs = (size_t)(end_pair.first - uniq_keys);
			BLST_ASSERT(runs <= (size_t)(W + 1));

			memcpyAsync<B>(
				d_runs,
				&runs,
				sizeof(size_t),
				MemcpyKind::HostToDevice,
				stream);

			GroupZZ* out_row = d_groups.data_handle() + (size_t)s * (size_t)W;

			auto ec = scatter_group_sums_multi<B, GroupZZ, ValT>(
				span<const uint32_t>(uniq_keys, W + 1),
				span<const ValT>(agg_vals, W + 1),
				d_runs,
				W,
				1,
				d2Span<GroupZZ>(out_row, 1, W),
				stream);
			if (ec != DeviceError::cudaSuccess)
				return ec;
		}

		return getLastError<B>();
	}


	// we have D groups, each with W values. for each group, we want to compute
	// the integrated sum of the values in reverse order.
	// 
	// High->low running sum over buckets(xyzz form) : 
	// computes sum_d d S_d from pre-summed buckets and writes 
	// the window's contribution as Jacobian.

	// we have D groups, each with W values. for each group, we want to compute
	// the integrated sum of the values in reverse order.
	// 
	// High->low running sum over buckets(xyzz form) : 
	// computes ?? d??S_d from pre-summed buckets and writes 
	// the window???s contribution as Jacobian.
	template<class Group, typename GroupZZ>
	BLST_HD BLST_INLINE void integrate_buckets_xyzz_kernel(
		int tid,
		const GroupZZ* BLST_restrict d_groups_xyzz,  // [D*W]
		const WindowExtract* BLST_restrict d_windows,
		int M,
		int D,
		int W,
		Group* BLST_restrict d_window_sums       // [M*D] Jacobian out
	)
	{
		if (tid >= D * M)
			return;
		auto d = tid / M;
		//auto m = tid % M;

		auto& window = d_windows[d];

		const GroupZZ* buckets = d_groups_xyzz + (size_t)tid * (size_t)W;

		Group out;
		{
			GroupZZ ret, acc;
			auto wbits = window.encode_size - 1;
			size_t n = (size_t)1 << wbits;

			/* Calculate sum of x[i-1]*i for i=1 through 1<<|wbits|. */
			vec_copy(&acc, &buckets[--n], sizeof(acc));
			vec_copy(&ret, &buckets[n], sizeof(ret));

			// computes Σ d·S_d
			while (n--) {
				GroupZZ::dadd(acc, acc, buckets[n]);
				GroupZZ::dadd(ret, ret, acc);
			}
			msm_to_jacobian(out, ret);
		}
		d_window_sums[tid] = out;
	}

#ifdef BLST_CUDA
	template<class Group, typename GroupZZ>
	BLST_GLOBAL void integrate_buckets_xyzz_kernel(
		const GroupZZ* BLST_restrict d_groups_xyzz,  // [D*W]
		const WindowExtract* BLST_restrict d_windows,
		int M,
		int D,
		int W,
		Group* BLST_restrict d_window_sums       // [M*D] Jacobian out
	)
	{
		int tid = blockDim.x * blockIdx.x + threadIdx.x;
		integrate_buckets_xyzz_kernel<Group, GroupZZ>(
			tid,
			d_groups_xyzz,
			d_windows,
			M,
			D,
			W,
			d_window_sums);
	}
#endif // CUDA


	template<Backend B,
		class Group, typename GroupZZ>
	DeviceError integrate_buckets_xyzz(
		d3Span<GroupZZ> d_groups_xyzz,  // [D*M*W]
		span<WindowExtract> d_windows, // [D]
		int M,
		int D,
		int W,
		d2Span<Group> d_window_sums,      // [D*M]
		Stream stream = 0)
	{
		BLST_ASSERT(d_groups_xyzz.extent(0) == (size_t)D);
		BLST_ASSERT(d_groups_xyzz.extent(1) == (size_t)M);
		BLST_ASSERT(d_groups_xyzz.extent(2) == (size_t)W);
		BLST_ASSERT(d_window_sums.extent(0) == (size_t)D);
		BLST_ASSERT(d_window_sums.extent(1) == (size_t)M);
		BLST_ASSERT(d_windows.size() == (size_t)D);
		auto total = M * D;
		if constexpr (B == Backend::CPU)
		{
			for (u64 i = 0; i < total; ++i)
			{
				integrate_buckets_xyzz_kernel<Group, GroupZZ>(
					i,
					d_groups_xyzz.data_handle(),
					d_windows.data(),
					M,
					D,
					W,
					d_window_sums.data_handle());
			}
			return DeviceError::cudaSuccess;
		}
		else
		{
#if defined(BLST_CUDA) && defined(BLST_CUDA_MSM)

			static_assert(B == Backend::CUDA);
			int blk = 256;
			int grd = (total + blk - 1) / blk;
			integrate_buckets_xyzz_kernel<Group, GroupZZ> BLST_LAUNCH(grd, blk, 0, stream) (
				d_groups_xyzz.data_handle(),
				d_windows.data(),
				M,
				D,
				W,
				d_window_sums.data_handle());
			return cudaGetLastError();
#endif // BLST_CUDA
		}
	}

	template<class Group>
	BLST_HD BLST_INLINE void fold_windows_kernel(
		int m,
		const Group* BLST_restrict d_window_sums, // [D*M], index = m*D + d
		int M,
		int D,
		int wbits,
		Group* BLST_restrict d_outputs            // [M]
	)
	{
		if (m >= M)
			return;

		Group acc;
		vec_zero(&acc, sizeof(acc));

		for (int d = D - 1; d >= 0; --d)
		{
			for (int j = 0; j < wbits; ++j)
				Group::doubling(acc, acc);

			const Group& win = d_window_sums[d * M + m];
			Group tmp;
			Group::dadd(tmp, acc, win);
			acc = tmp;
		}

		d_outputs[m] = acc;
	}

#ifdef BLST_CUDA
	template<typename Group>
	BLST_GLOBAL void fold_windows_kernel(
		const Group* BLST_restrict d_window_sums, // [M*D]
		int M,
		int D,
		int wbits,
		Group* BLST_restrict d_outputs            // [M]
	)
	{
		int m = blockDim.x * blockIdx.x + threadIdx.x;
		fold_windows_kernel<Group>(
			m,
			d_window_sums,
			M,
			D,
			wbits,
			d_outputs);
	}
#endif // BLST_CUDA

	template<Backend B, class Group>
	DeviceError fold_windows(
		d2Span<const Group> d_window_sums, // [D][M]
		int M,
		int D,
		int wbits,
		span<Group> d_outputs,           // [M]
		Stream stream = 0)
	{
		BLST_ASSERT(d_window_sums.extent(0) == (size_t)D);
		BLST_ASSERT(d_window_sums.extent(1) == (size_t)M);
		BLST_ASSERT(d_outputs.size() == (size_t)M);

		if constexpr (B == Backend::CPU)
		{
			for (u64 m = 0; m < (u64)M; m++)
			{
				fold_windows_kernel<Group>(
					(int)m,
					d_window_sums.data_handle(),
					M,
					D,
					wbits,
					d_outputs.data());
			}
			return DeviceError::cudaSuccess;
		}
		else
		{
#if defined(BLST_CUDA) && defined(BLST_CUDA_MSM)

			static_assert(B == Backend::CUDA);

			int blk = 256;
			int grd = (M + blk - 1) / blk;
			fold_windows_kernel<Group> BLST_LAUNCH(grd, blk, 0, stream) (
				d_window_sums.data_handle(),
				M,
				D,
				wbits,
				d_outputs.data());
			return cudaGetLastError();

#endif // BLST_CUDA
		}
	}

	// ===== Multi MSM pipeline: Booth, batched ==============================
	// Performs M simultaneous MSMs of N points each.
	// 
	// d_outputs[i] = sum_j d_points[i*N + j] * d_scalars[i*N + j]
	// 
	// Uses Booth encoding with wbits window size.
	// @d_workspace: temp buffer, or nullptr to query size
	// @workspace_bytes: in: size of d_workspace; out: required size if d_workspace is nullptr or too small
	// @d_points: size M*N array of points. The i'th MSM uses points [i*N .. (i+1)*N-1]
	// @d_scalars: size M*N scalars, each nbytes long, flat array. The i'th MSM uses scalars [i*N*nbytes .. (i+1)*N*nbytes-1]
	// @nbytes: number of bytes per scalar
	// @nbits: number of significant bits per scalar
	// @N: number of points per MSM
	// @M: number of MSMs
	// @wbits: window size in bits (Booth)
	// @d_outputs: size M array of output points
	// @stream: cuda stream
	template<
		Backend B,
		typename Group,
		typename GroupAffine = Group,
		typename GroupZZ = Group>
	DeviceError msm_booth_batched_pipeline(
		void* d_workspace,                          // nullptr => size query
		size_t& workspace_bytes,                    // in/out
		d2Span<const GroupAffine> d_points,  // [M][N]
		d3Span<const uint8_t> d_scalars,            // [M][N][nbytes],
		int nbytes,
		int nbits,
		int N,                                      // points per MSM
		int M,                                      // number of MSMs
		int wbits,                                  // window bits (Booth)
		/* final outputs: */
		span<Group> d_outputs,                    // [M]
		Stream stream = 0,
		bool debug = false,
		bool print = false)
	{
		BLST_ASSERT(d_points.extent(0) == (size_t)M);
		BLST_ASSERT(d_points.extent(1) == (size_t)N);
		BLST_ASSERT(d_scalars.extent(0) == (size_t)M);
		BLST_ASSERT(d_scalars.extent(1) == (size_t)N);
		BLST_ASSERT(d_scalars.extent(2) == (size_t)nbytes);
		BLST_ASSERT(d_outputs.size() == (size_t)M);


		if (N <= 0 || M <= 0 || wbits <= 0)
			return DeviceError::cudaErrorInvalidValue;

		// number of windows. If nbits % wbits == 0, we have an extra window for the overlap bits.
		int D = 1 + nbits / wbits;
		// Booth bucket count (no zero) {-W..-1,1..W}, 2^wbits - 1 values
		int W = 1 << (wbits - 1);

		// ---- Workspace sizing for downstream steps ----
		//size_t needed = 0;
		uint8_t* base = static_cast<uint8_t*>(d_workspace);
		uint8_t* ws = base;

		// buffer a will hold:
		// 1)
		//   - extracted digits:      u8[D][M][N]
		//   - extracted index_sign:  IdxSign[D][M][N]
		// 2) 
		//   - group summations:   GroupZZ[D][M][W]
		span<u8> buffer_a(ws, std::max<size_t>(
			(sizeof(u8) + sizeof(IdxSign)) * D * M * N + alignof(IdxSign),
			sizeof(GroupZZ) * D * M * W + alignof(GroupZZ)));
		ws += buffer_a.size();

		// buffer b will hold:
		// 1)
		//   - sorted digits:      u8[D][M][N]
		//   - sorted index_sign:  IdxSign[D][M][N]
		// 2)
		//   - window summations: Group[D][M]
		span<u8> buffer_b(ws, std::max<size_t>(
			(sizeof(u8) + sizeof(IdxSign)) * D * M * N + alignof(IdxSign),
			sizeof(Group) * D * M + alignof(Group)));
		ws += buffer_b.size();

		// Per-digit segmented sort workspace:
		size_t sort_ws = 0;
		{
			// Assume you invoke sort per MSM or per packed layout; here we size for one (D,N) pack.
			DeviceError ec = sort_digits_segmented<B>(
				nullptr,
				sort_ws,
				(u8*)nullptr,
				(u8*)nullptr,
				(IdxSign*) nullptr,
				(IdxSign*) nullptr,
				N,
				D * M,
				wbits,
				stream);
			if (ec != DeviceError::cudaSuccess)
				return ec;
		}

		// Per-digit group reduce workspace:
		size_t reduce_ws = 0;
		{
			// Size for one call of reduce_groups
			//Group* dummy_groups = nullptr;
			DeviceError ec = reduce_groups<B, GroupZZ>(
				nullptr,
				reduce_ws,
				{}, // keys
				(Defaulted<GroupZZ>*)nullptr,
				{}, // groups
				D * M,
				N,
				W,
				stream);
			if (ec != DeviceError::cudaSuccess && ec != DeviceError::cudaErrorInvalidValue)
				return ec;
		}
		span<u8> buffer_c(ws, std::max<size_t>(sort_ws, reduce_ws));
		ws += buffer_c.size();

		span<u8> buffer_d(ws, sizeof(WindowExtract) * D + alignof(WindowExtract));
		ws += buffer_d.size();

		size_t needed = static_cast<size_t>(ws - base);
		if (d_workspace == nullptr)
		{
			workspace_bytes = needed;
			return DeviceError::cudaSuccess;
		}
		if (workspace_bytes < needed)
		{
			workspace_bytes = needed;
			return DeviceError::cudaErrorInvalidValue;
		}

		// ---- Carve workspace ----

		// buffer a, shared between digit extraction and group summation
		auto buffer_a_iter = buffer_a.data();
		buffer_a_iter = align_up_ptr(buffer_a_iter, alignof(u8));
		d3Span<u8> d_digit(reinterpret_cast<u8*>(buffer_a_iter), D, M, N);    // [D][M][N]
		buffer_a_iter += sizeof(u8) * d_digit.size();

		buffer_a_iter = align_up_ptr(buffer_a_iter, alignof(IdxSign));
		d3Span<IdxSign> d_index_signs(reinterpret_cast<IdxSign*>(buffer_a_iter), D, M, N);    // [D][M][N]
		buffer_a_iter += sizeof(IdxSign) * d_index_signs.size();
		assert(buffer_a_iter <= buffer_a.data() + buffer_a.size());

		// reset buffer a iterator for group summation
		buffer_a_iter = buffer_a.data();
		buffer_a_iter = align_up_ptr(buffer_a_iter, alignof(GroupZZ));
		d3Span<GroupZZ> d_groups(reinterpret_cast<GroupZZ*>(buffer_a_iter), D, M, W); // [D][M][W]
		buffer_a_iter += d_groups.size() * sizeof(GroupZZ);
		assert(buffer_a_iter <= buffer_a.data() + buffer_a.size());

		// buffer b, for sorted digits/index_signs
		auto buffer_b_iter = buffer_b.data();
		buffer_b_iter = align_up_ptr(buffer_b_iter, alignof(u8));
		d3Span<u8> d_sorted_digit(buffer_b_iter, D, M, N);    // [D][M][N]
		buffer_b_iter += sizeof(u8) * d_sorted_digit.size();

		buffer_b_iter = align_up_ptr(buffer_b_iter, alignof(IdxSign));
		d3Span<IdxSign> d_sorted_index_signs(reinterpret_cast<IdxSign*>(buffer_b_iter), D, M, N);    // [D][M][N]
		buffer_b_iter += sizeof(IdxSign) * d_sorted_index_signs.size();
		assert(buffer_b_iter <= buffer_b.data() + buffer_b.size());

		// reset buffer b for window sums
		buffer_b_iter = buffer_b.data();
		buffer_b_iter = align_up_ptr(buffer_b_iter, alignof(Group));
		d2Span<Group> d_window_sums(reinterpret_cast<Group*>(buffer_b_iter), D, M); // [D][M]
		buffer_b_iter += d_window_sums.size() * sizeof(Group);
		assert(buffer_b_iter <= buffer_b.data() + buffer_b.size());

		// ---- Precompute window extract parameters ----
		std::vector<WindowExtract> windows(D);
		for (size_t i = 0; i < D; i++)
			windows[i] = WindowExtract(i, wbits, nbits);

		// copy windows to device
		auto buffer_d_iter = align_up_ptr(buffer_d.data(), alignof(WindowExtract));
		auto d_windows = span<WindowExtract>(
			reinterpret_cast<WindowExtract*>(buffer_d_iter),
			D);
		auto ec = memcpyAsync<B>(
			d_windows.data(),
			windows.data(),
			d_windows.size_bytes(),
			MemcpyKind::HostToDevice,
			stream);
		if (ec != DeviceError::cudaSuccess)
			return ec;

		// ---- Pipeline ----
		// 

		// helper to reinterpret d3Span<u8> as d2Span<u8>, combining the first two dims.
		auto d2Outter = []<typename T>(d3Span<T> span) -> d2Span<T>
		{
			return d2Span<T>(
				span.data_handle(),
				span.extent(0) * span.extent(1),
				span.extent(2));
		};


		// helper to reinterpret d3Span<u8> as d2Span<u8>, combining the second two dims.
		auto d2Inner = []<typename T>(d3Span<T> span) -> d2Span<T>
		{
			return d2Span<T>(
				span.data_handle(),
				span.extent(0),
				span.extent(1) * span.extent(2));
		};


		// we have M*N scalars, each nbytes long.
		// for each we need to extract D digits of wbits each.
		// The outer most dimention is the digit index (0..D-1).
		// the middle dimension is the MSM index (0..M-1).
		// the inner most dimension is the scalar index (0..N-1).
		ec = extract_digits_flat<B>(
			d2Outter(d_scalars),       // [M*N][nbytes]
			(M * N),         // total scalars
			d_windows,
			nbytes,
			nbits,
			wbits,
			D,
			d2Inner(d_digit),        // [D][M*N]  digit
			d2Inner(d_index_signs),  // [D][M*N]  index|sign
			stream
		);
		if (ec != DeviceError::cudaSuccess)
			return ec;

		//if (print)
		//    print_extract_digit<Field>(d_scalars, nbytes, nbits, wbits, N, M, D, d_digit, d_index_signs, stream);

		// for each of the D*M digit slices of N elements,
		// sort the digits and store the result in d_index_2d.
		// This will contain the gather
		ec = sort_digits_segmented<B>(
			buffer_c.data(),
			sort_ws,         // in/out; CUB will honor larger buffers
			d_digit.data_handle(),            // [D*N] per MSM or [D*M*N] packed
			d_sorted_digit.data_handle(),    // [D*N] per MSM or [D*M*N] packed
			d_index_signs.data_handle(),   // D*M vectors, each of length N. contains src idx
			d_sorted_index_signs.data_handle(),   // D*M vectors, each of length N. contains src idx
			N,
			D * M,
			wbits,
			stream);
		if (ec != DeviceError::cudaSuccess)
			return ec;

		//if (debug)
		//{
		//	// check that all points are reduced.
		//	size_t total_points = static_cast<size_t>(M) * static_cast<size_t>(N);
		//	auto redc_ok = all_of(
		//		stream,
		//		counting_iterator<size_t>(0),
		//		counting_iterator<size_t>(0) + total_points,
		//		[=] BLST_HD(size_t idx)
		//	{
		//		return
		//			Field::is_redc(d_points[idx].X) &&
		//			Field::is_redc(d_points[idx].Y);
		//	});
		//	if (!redc_ok)
		//	{
		//		printf("msm_booth_batched_pipeline: found non-reduced points\n");
		//		return DeviceError::cudaErrorInvalidValue;
		//	}
		//}

		//if (print)
		//    print_sort_digit_indexsigns<Field>(N, M, D, wbits, d_digit, d_index_signs, d_sorted_digit, d_sorted_index_signs, stream);

		// gather the points into Pointxyzz, applying sign.
		auto gather_iter = transform_iterator(
			d_sorted_index_signs.data_handle(),
			PackedToPointXYZZ<GroupAffine, Defaulted<GroupZZ>>{ d_points.data_handle(), M* N }
		);

		// now reduce each of the D*M digit slices of N points to W buckets.
		ec = reduce_groups<B, GroupZZ>(
			buffer_c.data(),
			reduce_ws,       // in/out; CUB will honor larger buffers
			d2Outter(d_sorted_digit),      // [D*M][N] 
			gather_iter,
			d2Outter(d_groups), // [D*M][W] (carve from workspace above)
			D * M,
			N,
			W,
			stream);
		if (ec != DeviceError::cudaSuccess)
			return ec;

		if (print)
			print_reduce_groups<B, Group, GroupAffine, GroupZZ>(N, M, D, W, d_sorted_digit, d_sorted_index_signs, d_groups, d_points, stream);

		ec = integrate_buckets_xyzz<B, Group, GroupZZ>(
			d_groups,
			d_windows,
			M,
			D,
			W,
			d_window_sums,
			stream);
		if (ec != DeviceError::cudaSuccess)
			return ec;

		if (print)
			print_integrate_buckets_xyzz<B, Group, GroupZZ>(M, D, W, wbits, d_groups, d_window_sums, stream);

		// 4) fold D window sums per MSM with wbits doublings in between
		ec = fold_windows<B, Group>(
			d_window_sums,
			M,
			D,
			wbits,
			d_outputs,
			stream);
		if (ec != DeviceError::cudaSuccess)
			return ec;

		if (print)
			print_fold_windows<B, Group>(M, D, wbits, d_window_sums, d_outputs, stream);

		return DeviceError::cudaSuccess;
	}


}
