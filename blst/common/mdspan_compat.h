#pragma once

#include <cstddef>
#include <type_traits>

namespace blst::mdspan_compat
{
    template <class IndexT, std::size_t Rank>
    struct dextents
    {
        static_assert(std::is_integral_v<IndexT>);
        static constexpr std::size_t rank() { return Rank; }
    };

    template <class T, class Extents>
    class mdspan;

    template <class T, class IndexT>
    class mdspan<T, dextents<IndexT, 2>>
    {
    public:
        using element_type = T;
        using index_type = IndexT;

        mdspan() = default;

        mdspan(T* ptr, IndexT e0, IndexT e1)
            : ptr_(ptr)
            , e0_(e0)
            , e1_(e1)
        {
        }

        template <class U,
            class = std::enable_if_t<std::is_convertible_v<U*, T*>>>
        mdspan(const mdspan<U, dextents<IndexT, 2>>& other)
            : ptr_(other.data_handle())
            , e0_(other.extent(0))
            , e1_(other.extent(1))
        {
        }

        T* data_handle() const
        {
            return ptr_;
        }

        IndexT extent(std::size_t i) const
        {
            return (i == 0) ? e0_ : e1_;
        }

        std::size_t size() const
        {
            return static_cast<std::size_t>(e0_) * static_cast<std::size_t>(e1_);
        }

        T& operator()(IndexT i0, IndexT i1) const
        {
            return ptr_[static_cast<std::size_t>(i0) * static_cast<std::size_t>(e1_) +
                static_cast<std::size_t>(i1)];
        }

    private:
        T* ptr_ = nullptr;
        IndexT e0_ = 0;
        IndexT e1_ = 0;
    };

    template <class T, class IndexT>
    class mdspan<T, dextents<IndexT, 3>>
    {
    public:
        using element_type = T;
        using index_type = IndexT;

        mdspan() = default;

        mdspan(T* ptr, IndexT e0, IndexT e1, IndexT e2)
            : ptr_(ptr)
            , e0_(e0)
            , e1_(e1)
            , e2_(e2)
        {
        }

        template <class U,
            class = std::enable_if_t<std::is_convertible_v<U*, T*>>>
        mdspan(const mdspan<U, dextents<IndexT, 3>>& other)
            : ptr_(other.data_handle())
            , e0_(other.extent(0))
            , e1_(other.extent(1))
            , e2_(other.extent(2))
        {
        }

        T* data_handle() const
        {
            return ptr_;
        }

        IndexT extent(std::size_t i) const
        {
            if (i == 0)
            {
                return e0_;
            }
            if (i == 1)
            {
                return e1_;
            }
            return e2_;
        }

        std::size_t size() const
        {
            return static_cast<std::size_t>(e0_) *
                static_cast<std::size_t>(e1_) *
                static_cast<std::size_t>(e2_);
        }

        T& operator()(IndexT i0, IndexT i1, IndexT i2) const
        {
            std::size_t s1 = static_cast<std::size_t>(e1_);
            std::size_t s2 = static_cast<std::size_t>(e2_);
            std::size_t off =
                static_cast<std::size_t>(i0) * (s1 * s2) +
                static_cast<std::size_t>(i1) * s2 +
                static_cast<std::size_t>(i2);
            return ptr_[off];
        }

    private:
        T* ptr_ = nullptr;
        IndexT e0_ = 0;
        IndexT e1_ = 0;
        IndexT e2_ = 0;
    };

    template <class T>
    using d2_span = mdspan<T, dextents<std::size_t, 2>>;

    template <class T>
    using d3_span = mdspan<T, dextents<std::size_t, 3>>;
} // namespace blst::mdspan_compat
