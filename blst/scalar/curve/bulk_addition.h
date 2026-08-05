/*
 * Copyright Supranational LLC
 * Copyright Category Labs
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#include "blst/fields/Fp.h"
#include "blst/fields/Fp2.h"
#include "blst/scalar/curve/point.h"

#include <cstddef>
#include <new>
#include <type_traits>

namespace blst::scalar
{

#ifndef BLST_SCALAR_BULK_SCRATCH_LIMIT
    /* Performance with 144 KiB scratch is within 1-3% of optimal. */
# define BLST_SCALAR_BULK_SCRATCH_LIMIT (144 * 1024)
#endif

    namespace detail
    {
        /*
         * Scratch for batched affine addition. `work` holds the product chain
         * and inverse denominator; it is not a Jacobian Z coordinate. Keeping
         * the affine value as a real subobject makes every mixed-add operand a
         * valid C++ object without changing the historic x|y|work layout.
         */
        template<typename Field>
        struct BatchAffineNode
        {
            AffinePoint<Field> point;
            Field work;
        };

        template<typename Field>
        inline constexpr bool batch_affine_node_layout =
            std::is_standard_layout_v<BatchAffineNode<Field>> &&
            std::is_trivially_default_constructible_v<BatchAffineNode<Field>> &&
            std::is_trivially_copyable_v<BatchAffineNode<Field>> &&
            std::is_trivially_destructible_v<BatchAffineNode<Field>> &&
            sizeof(BatchAffineNode<Field>) == sizeof(Point<Field>) &&
            alignof(BatchAffineNode<Field>) == alignof(Point<Field>) &&
            alignof(BatchAffineNode<Field>) <= alignof(limb_t) &&
            sizeof(BatchAffineNode<Field>) % sizeof(limb_t) == 0 &&
            offsetof(BatchAffineNode<Field>, point) == 0 &&
            offsetof(BatchAffineNode<Field>, work) == offsetof(Point<Field>, z);

        static_assert(batch_affine_node_layout<Fp>);
        static_assert(batch_affine_node_layout<Fp2>);

        /*
         * A C caller supplies fixed-table scratch as limb storage. Constructing
         * each node at its byte offset starts the C++ object's lifetime before
         * the gather writes it. Heap/vector callers use the same zero-overhead
         * view. All nodes used by points_accumulate have already been started.
         */
        template<typename Field>
        class BatchAffineWorkspace
        {
        public:
            using Node = BatchAffineNode<Field>;

            BLST_HD constexpr BatchAffineWorkspace() : storage_(nullptr), size_(0) {}
            BLST_HD constexpr BatchAffineWorkspace(void* storage, size_t size)
                : storage_(storage), size_(size) {}

            BLST_HD constexpr bool empty() const { return size_ == 0; }
            BLST_HD constexpr size_t size() const { return size_; }

            BLST_HD BLST_INLINE Node& start(size_t index) const
            {
                BLST_ASSERT(index < size_);
                byte* address = static_cast<byte*>(storage_) + index * sizeof(Node);
                BLST_ASSERT(reinterpret_cast<uintptr_t>(address) % alignof(Node) == 0);
                return *::new (static_cast<void*>(address)) Node;
            }

            BLST_HD BLST_INLINE Node& get(size_t index) const
            {
                BLST_ASSERT(index < size_);
                byte* address = static_cast<byte*>(storage_) + index * sizeof(Node);
                BLST_ASSERT(reinterpret_cast<uintptr_t>(address) % alignof(Node) == 0);
                return *std::launder(reinterpret_cast<Node*>(address));
            }

        private:
            void* storage_;
            size_t size_;
        };

        struct BulkScratchAllocator
        {
            static void* allocate(size_t bytes) noexcept
            {
                return ::operator new(bytes, std::nothrow);
            }

            static void deallocate(void* storage) noexcept
            {
                ::operator delete(storage);
            }
        };

        /*
         * This variable-time implementation amortizes field inversion across
         * affine additions using Montgomery's batch-inversion method. It
         * handles equal, opposite, and infinity inputs through the complete
         * mixed-addition cases below. Inputs and schedules must be public.
         */
        template<typename Field>
        BLST_HD BLST_INLINE void point_head(BatchAffineNode<Field>& A,
                                            BatchAffineNode<Field>& B,
                                            const Field* mul_acc)
        {
            limb_t inf = vec_is_zero(&A.point, sizeof(A.point)) |
                         vec_is_zero(&B.point, sizeof(B.point));
            const Field zero = Field::zero();

            Field::sub(B.work, B.point.x, A.point.x);       /* X2-X1 */
            Field::add(B.point.x, B.point.x, A.point.x);   /* X2+X1 */
            Field::add(A.work, B.point.y, A.point.y);      /* Y2+Y1 */
            Field::sub(B.point.y, B.point.y, A.point.y);   /* Y2-Y1 */
            if (vec_is_zero(&B.work, sizeof(B.work))) {    /* X2==X1 */
                inf = vec_is_zero(&A.work, sizeof(A.work));
                vec_select(&B.point.x, &A.work, &B.point.x, sizeof(B.point.x), inf);
                Field::sqr(B.point.y, A.point.x);
                Field::mul_by_3(B.point.y, B.point.y);     /* 3*X1^2 */
                vec_copy(&B.work, &A.work, sizeof(B.work));/* 2*Y1 */
            }                                             /* B.y numerator */
                                                          /* B.work denominator */
            vec_select(&A.point.x, &B.point.x, &A.point.x, sizeof(A.point.x), inf);
            vec_select(&A.point.y, &A.work, &A.point.y, sizeof(A.point.y), inf);
            vec_select(&A.work, &Field::one(), &B.work, sizeof(A.work), inf);
            vec_select(&B.work, &zero, &B.work, sizeof(B.work), inf);
            if (mul_acc != nullptr)
                Field::mul(A.work, A.work, *mul_acc);
        }

        template<typename Field>
        BLST_HD BLST_INLINE void point_tail(BatchAffineNode<Field>& D,
                                            BatchAffineNode<Field>& A,
                                            BatchAffineNode<Field>& B,
                                            Field& lambda)
        {
            Field lambda_squared;
            limb_t inf = vec_is_zero(&B.work, sizeof(B.work));

            Field::mul(lambda, lambda, B.point.y);
            Field::sqr(lambda_squared, lambda);
            Field::sub(D.point.x, lambda_squared, B.point.x);

            Field::sub(D.point.y, A.point.x, D.point.x);
            Field::mul(D.point.y, D.point.y, lambda);
            Field::sub(D.point.y, D.point.y, A.point.y);

            vec_select(&D.point.x, &A.point.x, &D.point.x,
                       2 * sizeof(D.point.x), inf);
            vec_select(&B.work, &Field::one(), &B.work, sizeof(B.work), inf);
        }

        /* Add nodes[0..n) to sum. Every node must be a live object. */
        template<typename Field>
        BLST_HD BLST_INLINE void points_accumulate(
            Point<Field>& sum,
            const BatchAffineWorkspace<Field>& workspace,
            size_t n)
        {
            size_t base = 0;

            while (n >= 16) {
                if (n & 1) {
                    Point<Field>::dadd(sum, sum, workspace.get(base).point);
                    ++base;
                }
                n /= 2;
                for (size_t i = 0; i < n; ++i) {
                    auto& A = workspace.get(base + 2 * i);
                    auto& B = workspace.get(base + 2 * i + 1);
                    const Field* mul_acc = i == 0
                        ? nullptr
                        : &workspace.get(base + 2 * i - 2).work;
                    point_head<Field>(A, B, mul_acc);
                }
                auto& last_work = workspace.get(base + 2 * n - 2).work;
                Field::reciprocal(last_work, last_work);

                for (size_t i = n - 1; i != 0; --i) {
                    auto& A = workspace.get(base + 2 * i);
                    auto& B = workspace.get(base + 2 * i + 1);
                    auto& previous_A = workspace.get(base + 2 * i - 2);
                    auto& destination = workspace.get(base + n + i);
                    Field::mul(previous_A.work, A.work, previous_A.work);
                    point_tail<Field>(destination, A, B, previous_A.work);
                    Field::mul(previous_A.work, A.work, B.work);
                }
                point_tail<Field>(workspace.get(base + n),
                                  workspace.get(base),
                                  workspace.get(base + 1),
                                  workspace.get(base).work);
                base += n;
            }
            while (n--)
                Point<Field>::dadd(sum, sum, workspace.get(base++).point);
        }

        template<typename Field>
        BLST_INLINE void points_add_serial(Point<Field>& sum,
                                           const AffinePoint<Field>* const points[],
                                           size_t npoints)
        {
            const AffinePoint<Field>* point = nullptr;
            while (npoints--) {
                point = *points ? *points++ : point + 1;
                Point<Field>::dadd(sum, sum, *point);
            }
        }
    }

    /**
     * @brief Add public subgroup affine points using bounded batch inversion.
     * @param sum Receives the complete sum; exact aliasing with input storage is
     * unsupported.
     * @param points Explicit pointers or `{contiguous_base, nullptr}` compact form.
     * @param npoints Number of affine points represented by @p points.
     * @post Invalid/empty input produces infinity. Allocation failure falls back
     * to complete serial addition and never publishes a partial result.
     * @note Variable-time. The one bounded allocation is capped by
     * `BLST_SCALAR_BULK_SCRATCH_LIMIT` and is skipped for fewer than 16 points.
     */
    template<typename Field, typename Allocator = detail::BulkScratchAllocator>
    BLST_INLINE void points_add(Point<Field>& sum,
                                const AffinePoint<Field>* const points[],
                                size_t npoints)
    {
        vec_zero(&sum, sizeof(sum));
        if (npoints == 0 || points == nullptr || points[0] == nullptr)
            return;

        if (npoints < 16) {
            detail::points_add_serial(sum, points, npoints);
            return;
        }

        using Node = detail::BatchAffineNode<Field>;
        constexpr size_t max_nodes = BLST_SCALAR_BULK_SCRATCH_LIMIT / sizeof(Node);
        static_assert(max_nodes >= 16,
                      "BLST_SCALAR_BULK_SCRATCH_LIMIT must hold at least 16 nodes");
        const size_t capacity = npoints < max_nodes ? npoints : max_nodes;
        void* storage = Allocator::allocate(capacity * sizeof(Node));
        if (storage == nullptr) {
            detail::points_add_serial(sum, points, npoints);
            return;
        }

        detail::BatchAffineWorkspace<Field> workspace(storage, capacity);
        const AffinePoint<Field>* point = nullptr;
        while (npoints != 0) {
            const size_t count = npoints < capacity ? npoints : capacity;
            for (size_t i = 0; i < count; ++i) {
                point = *points ? *points++ : point + 1;
                auto& node = workspace.start(i);
                vec_copy(&node.point, point, sizeof(node.point));
            }
            detail::points_accumulate(sum, workspace, count);
            npoints -= count;
        }
        Allocator::deallocate(storage);
    }
}
