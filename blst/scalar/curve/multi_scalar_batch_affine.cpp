/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/scalar/curve/multi_scalar.h"

#include <cstdint>
#include <limits>
#include <new>
#include <type_traits>

namespace blst::scalar
{
	namespace
	{
		template<typename Field>
		struct BatchAffineBucket
		{
			Field x;
			Field y;
			std::uint32_t generation;
			std::uint32_t queued_epoch;
		};

		template<typename Field>
		struct BatchProjectiveBucket
		{
			Point<Field> point;
			std::uint32_t generation;
		};

		struct BatchSchedulePoint
		{
			std::uint32_t point_index;
			std::uint32_t bucket_index;
			std::uint32_t subtract;
		};

		template<typename T>
		class ScratchArray
		{
			byte* storage_;
		#ifndef NDEBUG
			size_t size_;
		#endif

			BLST_INLINE byte* address(size_t index) const
			{
		#ifndef NDEBUG
				BLST_ASSERT(index < size_);
		#endif
				return storage_ + index * sizeof(T);
			}

		public:
			ScratchArray() : storage_(nullptr)
		#ifndef NDEBUG
				, size_(0)
		#endif
			{}
			ScratchArray(byte* storage, size_t size)
				: storage_(storage)
		#ifndef NDEBUG
				, size_(size)
		#endif
			{
				static_assert(std::is_trivially_default_constructible_v<T>);
				static_assert(std::is_trivially_destructible_v<T>);
				(void)size;
				BLST_ASSERT(
					reinterpret_cast<uintptr_t>(storage_) % alignof(T) == 0);
			}

			void start_all(size_t count) const
			{
				BLST_ASSERT(count == size_);
				for (size_t i = 0; i < count; ++i)
					::new (static_cast<void*>(storage_ + i * sizeof(T))) T;
			}

			BLST_INLINE T& operator[](size_t index)
			{
				return *std::launder(reinterpret_cast<T*>(address(index)));
			}

			BLST_INLINE const T& operator[](size_t index) const
			{
				return *std::launder(reinterpret_cast<const T*>(address(index)));
			}
		};

		template<typename Field, size_t BatchSize>
		struct BatchAffineWorkspace
		{

			ScratchArray<BatchAffineBucket<Field>> affine_buckets;
			ScratchArray<BatchProjectiveBucket<Field>> collision_buckets;
			ScratchArray<BatchSchedulePoint> schedule;
			ScratchArray<Field> numerators;
			ScratchArray<Field> denominators;
			size_t bucket_count;
			size_t schedule_size;
			std::uint32_t window_generation;
			std::uint32_t queue_epoch;

			void start_lifetimes() const
			{
				affine_buckets.start_all(bucket_count);
				collision_buckets.start_all(bucket_count);
				schedule.start_all(BatchSize);
				numerators.start_all(BatchSize);
				denominators.start_all(BatchSize);
			}
		};

		constexpr bool checked_add(size_t left, size_t right, size_t& out)
		{
			if (right > std::numeric_limits<size_t>::max() - left)
				return false;
			out = left + right;
			return true;
		}

		constexpr bool checked_align(size_t value, size_t alignment, size_t& out)
		{
			BLST_ASSERT(alignment != 0 && (alignment & (alignment - 1)) == 0);
			size_t with_padding;
			if (!checked_add(value, alignment - 1, with_padding))
				return false;
			out = with_padding & ~(alignment - 1);
			return true;
		}

		template<typename T>
		constexpr bool add_array(size_t count, size_t& offset)
		{
			size_t aligned;
			size_t bytes;
			if (!checked_align(offset, alignof(T), aligned) ||
				!detail::checked_mul_size(count, sizeof(T), bytes) ||
				!checked_add(aligned, bytes, offset))
				return false;
			return true;
		}

		template<typename Field, size_t BatchSize>
		size_t batch_affine_scratch_size(size_t window_size)
		{
			static_assert(BatchSize > 0);
			static_assert(alignof(BatchAffineBucket<Field>) <= alignof(limb_t));
			static_assert(alignof(BatchProjectiveBucket<Field>) <= alignof(limb_t));
			static_assert(alignof(BatchSchedulePoint) <= alignof(limb_t));
			static_assert(alignof(Field) <= alignof(limb_t));

			if (!detail::valid_pippenger_window(window_size))
				return 0;

			const size_t bucket_count = size_t{1} << (window_size - 1);
			size_t offset = 0;
			if (!add_array<BatchAffineBucket<Field>>(bucket_count, offset) ||
				!add_array<BatchProjectiveBucket<Field>>(bucket_count, offset) ||
				!add_array<BatchSchedulePoint>(BatchSize, offset) ||
				!add_array<Field>(BatchSize, offset) ||
				!add_array<Field>(BatchSize, offset) ||
				!checked_align(offset, alignof(limb_t), offset))
				return 0;
			return offset;
		}

		struct ScratchCursor
		{
			byte* base;
			size_t offset;

			template<typename T>
			auto take(size_t count)
			{
				static_assert((alignof(T) & (alignof(T) - 1)) == 0);
				BLST_ASSERT(
					offset <= std::numeric_limits<size_t>::max() - (alignof(T) - 1));
				const size_t aligned =
					(offset + alignof(T) - 1) & ~(alignof(T) - 1);
				BLST_ASSERT(
					count <= (std::numeric_limits<size_t>::max() - aligned) /
						sizeof(T));
				offset = aligned + count * sizeof(T);
				return ScratchArray<T>{
					base + aligned, count};
			}
		};

		template<typename Field, size_t BatchSize>
		BatchAffineWorkspace<Field, BatchSize> make_workspace(
			limb_t scratch[],
			size_t bucket_count)
		{
			ScratchCursor cursor{
				reinterpret_cast<byte*>(scratch),
				0,
			};
			BatchAffineWorkspace<Field, BatchSize> workspace{
				cursor.take<BatchAffineBucket<Field>>(bucket_count),
				cursor.take<BatchProjectiveBucket<Field>>(bucket_count),
				cursor.take<BatchSchedulePoint>(BatchSize),
				cursor.take<Field>(BatchSize),
				cursor.take<Field>(BatchSize),
				bucket_count,
				0,
				0,
				0,
			};
			workspace.start_lifetimes();
			return workspace;
		}

		template<typename Field>
		BLST_INLINE bool reduced_equal(const Field& left, const Field& right)
		{
			if constexpr (std::is_same_v<Field, Fp>) {
				BLST_ASSERT(Fp::is_redc(left));
				BLST_ASSERT(Fp::is_redc(right));
			}
			else {
				static_assert(std::is_same_v<Field, Fp2>);
				BLST_ASSERT(Fp::is_redc(left[0]));
				BLST_ASSERT(Fp::is_redc(left[1]));
				BLST_ASSERT(Fp::is_redc(right[0]));
				BLST_ASSERT(Fp::is_redc(right[1]));
			}
			return vec_is_equal(&left, &right, sizeof(Field)) != 0;
		}

		template<typename Field>
		BLST_INLINE bool field_is_zero(const Field& value)
		{
			return vec_is_zero(&value, sizeof(value)) != 0;
		}

		template<typename Field, size_t BatchSize>
		void advance_window_generation(
			BatchAffineWorkspace<Field, BatchSize>& workspace)
		{
			if (++workspace.window_generation != 0)
				return;

			for (size_t i = 0; i < workspace.bucket_count; ++i) {
				workspace.affine_buckets[i].generation = 0;
				workspace.collision_buckets[i].generation = 0;
			}
			workspace.window_generation = 1;
		}

		template<typename Field, size_t BatchSize>
		void advance_queue_epoch(
			BatchAffineWorkspace<Field, BatchSize>& workspace)
		{
			if (++workspace.queue_epoch != 0)
				return;

			for (size_t i = 0; i < workspace.bucket_count; ++i)
				workspace.affine_buckets[i].queued_epoch = 0;
			workspace.queue_epoch = 1;
		}

		template<typename Field, size_t BatchSize>
		void batch_add(
			BatchAffineWorkspace<Field, BatchSize>& workspace,
			const AffinePoint<Field> points[])
		{
			Field accumulator = Field::one();

			for (size_t i = 0; i < workspace.schedule_size; ++i) {
				const BatchSchedulePoint item = workspace.schedule[i];
				BatchAffineBucket<Field>& bucket =
					workspace.affine_buckets[item.bucket_index];
				const AffinePoint<Field>& point = points[item.point_index];

				Field signed_y;
				Field::cneg(signed_y, point.y, item.subtract != 0);
				if (reduced_equal(bucket.x, point.x)) {
					if (!reduced_equal(bucket.y, signed_y)) {
						bucket.generation = 0;
						continue;
					}

					Field x_squared;
					Field::sqr(x_squared, point.x);
					Field::mul_by_3(x_squared, x_squared);
					Field::add(
						workspace.denominators[i], bucket.y, bucket.y);
					if (field_is_zero(workspace.denominators[i])) {
						bucket.generation = 0;
						continue;
					}
					Field::mul(
						workspace.numerators[i], accumulator, x_squared);
					Field::mul(
						accumulator, accumulator,
						workspace.denominators[i]);
					continue;
				}

				Field numerator;
				Field::sub(workspace.denominators[i], bucket.x, point.x);
				Field::sub(numerator, bucket.y, signed_y);
				Field::mul(
					workspace.numerators[i], accumulator, numerator);
				Field::mul(
					accumulator, accumulator, workspace.denominators[i]);
			}

			Field::reciprocal(accumulator, accumulator);

			for (size_t i = workspace.schedule_size; i-- > 0;) {
				const BatchSchedulePoint item = workspace.schedule[i];
				BatchAffineBucket<Field>& bucket =
					workspace.affine_buckets[item.bucket_index];
				if (bucket.generation != workspace.window_generation)
					continue;
				const AffinePoint<Field>& point = points[item.point_index];

				Field lambda;
				Field::mul(lambda, accumulator, workspace.numerators[i]);
				Field::mul(
					accumulator, accumulator, workspace.denominators[i]);

				Field x;
				Field::sqr(x, lambda);
				Field x_sum;
				Field::add(x_sum, bucket.x, point.x);
				Field::sub(x, x, x_sum);

				Field y;
				Field::sub(y, point.x, x);
				Field::mul(y, lambda, y);
				Field signed_y;
				Field::cneg(signed_y, point.y, item.subtract != 0);
				Field::sub(y, y, signed_y);

				bucket.x = x;
				bucket.y = y;
			}
		}

		template<typename Field, size_t BatchSize>
		void flush_schedule(
			BatchAffineWorkspace<Field, BatchSize>& workspace,
			const AffinePoint<Field> points[])
		{
			if (workspace.schedule_size == 0)
				return;
			batch_add(workspace, points);
			workspace.schedule_size = 0;
			advance_queue_epoch(workspace);
		}

		template<typename Field, size_t BatchSize>
		void add_projective_collision(
			BatchAffineWorkspace<Field, BatchSize>& workspace,
			size_t bucket_index,
			const AffinePoint<Field>& point,
			bool subtract)
		{
			BatchProjectiveBucket<Field>& bucket =
				workspace.collision_buckets[bucket_index];
			AffinePoint<Field> signed_point = point;
			AffinePoint<Field>::cneg(signed_point, subtract);
			if (bucket.generation != workspace.window_generation) {
				Point<Field>::from_affine(bucket.point, signed_point);
				bucket.generation = workspace.window_generation;
				return;
			}

			Point<Field> result;
			Point<Field>::dadd(result, bucket.point, signed_point);
			bucket.point = result;
		}

		template<typename Field, size_t BatchSize>
		void add_to_bucket(
			BatchAffineWorkspace<Field, BatchSize>& workspace,
			const AffinePoint<Field> points[],
			size_t point_index,
			size_t bucket_index,
			bool subtract)
		{
			BatchAffineBucket<Field>& bucket =
				workspace.affine_buckets[bucket_index];
			if (bucket.generation != workspace.window_generation) {
				bucket.x = points[point_index].x;
				Field::cneg(bucket.y, points[point_index].y, subtract);
				bucket.generation = workspace.window_generation;
				return;
			}

			if (bucket.queued_epoch == workspace.queue_epoch) {
				add_projective_collision(
					workspace, bucket_index, points[point_index], subtract);
				return;
			}

			workspace.schedule[workspace.schedule_size++] =
				BatchSchedulePoint{
					static_cast<std::uint32_t>(point_index),
					static_cast<std::uint32_t>(bucket_index),
					subtract ? 1u : 0u,
				};
			bucket.queued_epoch = workspace.queue_epoch;
			if (workspace.schedule_size == BatchSize)
				flush_schedule(workspace, points);
		}

		template<typename Field, size_t BatchSize>
		Point<Field> reduce_buckets(
			const BatchAffineWorkspace<Field, BatchSize>& workspace,
			size_t bucket_count)
		{
			Pointxyzz<Field> running = {};
			Pointxyzz<Field> total = {};
			for (size_t bucket_index = bucket_count; bucket_index-- > 0;) {
				const BatchAffineBucket<Field>& affine =
					workspace.affine_buckets[bucket_index];
				const BatchProjectiveBucket<Field>& collision =
					workspace.collision_buckets[bucket_index];

				Pointxyzz<Field> bucket_sum = {};
				const bool has_affine =
					affine.generation == workspace.window_generation;
				const bool has_collision =
					collision.generation == workspace.window_generation;
				if (has_collision) {
					bucket_sum.x = collision.point.x;
					bucket_sum.y = collision.point.y;
					Field::sqr(bucket_sum.zz, collision.point.z);
					Field::mul(
						bucket_sum.zzz, bucket_sum.zz,
						collision.point.z);
				}
				if (has_affine) {
					const AffinePoint<Field> point{affine.x, affine.y};
					if (has_collision) {
						Pointxyzz<Field> result;
						Pointxyzz<Field>::dadd(
							result, bucket_sum, point, false);
						bucket_sum = result;
					}
					else {
						Pointxyzz<Field>::from_affine(bucket_sum, point);
					}
				}

				Pointxyzz<Field> result;
				Pointxyzz<Field>::dadd(result, running, bucket_sum);
				running = result;
				Pointxyzz<Field>::dadd(result, total, running);
				total = result;
			}

			Point<Field> out;
			Pointxyzz<Field>::to_Jacobian(out, total);
			return out;
		}
	}

	template<typename Field, size_t BatchSize>
	size_t msm_mult_batch_affine_scratch_sizeof(size_t window_size)
	{
		return batch_affine_scratch_size<Field, BatchSize>(window_size);
	}

	template<typename Field, size_t BatchSize>
	bool msm_mult_batch_affine(
		Point<Field>& ret,
		const AffinePoint<Field> points[],
		const byte scalars[],
		size_t npoints,
		size_t nbits,
		limb_t scratch[],
		size_t scratch_size,
		size_t window_size)
	{
		if (npoints == 0) {
			ret = {};
			return true;
		}
		if (points == nullptr || scalars == nullptr || scratch == nullptr ||
			npoints > std::numeric_limits<std::uint32_t>::max() ||
			nbits == 0 || nbits > 256)
			return false;

		const size_t required =
			batch_affine_scratch_size<Field, BatchSize>(window_size);
		if (required == 0 || scratch_size < required)
			return false;

		const size_t bucket_count = size_t{1} << (window_size - 1);
		vec_zero(scratch, required);
		auto workspace =
			make_workspace<Field, BatchSize>(scratch, bucket_count);
		BLST_ASSERT(workspace.bucket_count == bucket_count);
		for (size_t i = 0; i < npoints; ++i) {
			BLST_ASSERT(!AffinePoint<Field>::is_inf(points[i]));
			BLST_ASSERT(AffinePoint<Field>::in_subgroup(points[i]));
		}
		const size_t scalar_bytes = (nbits + 7) / 8;
		const size_t window_count = 1 + nbits / window_size;
		ret = {};

		for (size_t window_index = window_count; window_index-- > 0;) {
			advance_window_generation(workspace);
			advance_queue_epoch(workspace);
			workspace.schedule_size = 0;

			const WindowExtract extract(
				window_index, window_size, nbits);
			const size_t active_bucket_count =
				size_t{1} << (extract.encode_size - 1);
			const limb_t digit_mask =
				(limb_t{1} << extract.encode_size) - 1;
			for (size_t point_index = 0;
				point_index < npoints; ++point_index) {
				const limb_t encoding = extract.extract(
					scalars + point_index * scalar_bytes);
				limb_t digit = encoding & digit_mask;
				if (digit == 0)
					continue;
				const bool subtract =
					((encoding >> extract.encode_size) & 1) != 0;
				--digit;
				add_to_bucket(
					workspace, points, point_index,
					static_cast<size_t>(digit), subtract);
			}
			flush_schedule(workspace, points);

			const Point<Field> tile =
				reduce_buckets(workspace, active_bucket_count);
			Point<Field> combined;
			Point<Field>::dadd(combined, ret, tile);
			ret = combined;
			if (window_index != 0) {
				for (size_t bit = 0; bit < window_size; ++bit) {
					Point<Field> doubled;
					Point<Field>::doubling(doubled, ret);
					ret = doubled;
				}
			}
		}
		return true;
	}

#define BLST_INSTANTIATE_BATCH_AFFINE(Field, BatchSize) \
	template size_t msm_mult_batch_affine_scratch_sizeof<Field, BatchSize>(size_t); \
	template bool msm_mult_batch_affine<Field, BatchSize>( \
		Point<Field>&, const AffinePoint<Field>[], const byte[], size_t, \
		size_t, limb_t[], size_t, size_t)

	BLST_INSTANTIATE_BATCH_AFFINE(Fp, 128);
	BLST_INSTANTIATE_BATCH_AFFINE(Fp, 384);
	BLST_INSTANTIATE_BATCH_AFFINE(Fp2, 128);
	BLST_INSTANTIATE_BATCH_AFFINE(Fp2, 384);

#undef BLST_INSTANTIATE_BATCH_AFFINE
}
