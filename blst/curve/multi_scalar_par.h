#pragma once

#include "multi_scalar.h"

#ifdef BLST_CUDA
#include <cuda_runtime.h>

#ifndef __INTELLISENSE__
#include <thrust/iterator/iterator_facade.h>

#include <cub/cub.cuh>
#include <cub/device/device_reduce.cuh>
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


namespace blst
{

	enum class Backend {
		CPU
#if defined(BLST_CUDA)
		, CUDA
#endif
	};


#if defined(BLST_CUDA) && !defined(__INTELLISENSE__)

	static constexpr Backend DefaultBackend = Backend::CUDA;
	template<typename T>
	using d2Span = mdspan_compat::d2_span<T>;

	template<typename T>
	using d3Span = mdspan_compat::d3_span<T>;
#else

	static constexpr Backend DefaultBackend = Backend::CPU;
	template<typename T>
	using d2Span = mdspan_compat::d2_span<T>;
	template<typename T>
	using d3Span = mdspan_compat::d3_span<T>;
#endif


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
		Stream stream);


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
		Stream                     stream);

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
		bool print = false);

}
