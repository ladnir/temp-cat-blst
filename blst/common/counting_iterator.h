#include <iterator>
#include <cstdint>
#include <type_traits>


namespace blst
{

    template<class T>
    struct counting_iterator
    {
        static_assert(std::is_integral_v<T>);

        using value_type = T;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::random_access_iterator_tag;
        using reference = T;
        using pointer = void;

        T cur;

        counting_iterator() = default;
        explicit counting_iterator(T start) : cur(start) {}

        reference operator*() const
        {
            return cur;
        }

        counting_iterator& operator++()
        {
            cur += 1;
            return *this;
        }

        counting_iterator operator++(int)
        {
            counting_iterator tmp{ *this };
            ++(*this);
            return tmp;
        }

        counting_iterator& operator--()
        {
            cur -= 1;
            return *this;
        }

        counting_iterator operator--(int)
        {
            counting_iterator tmp{ *this };
            --(*this);
            return tmp;
        }

        counting_iterator& operator+=(difference_type n)
        {
            cur += static_cast<T>(n);
            return *this;
        }

        counting_iterator& operator-=(difference_type n)
        {
            cur -= static_cast<T>(n);
            return *this;
        }

        reference operator[](difference_type n) const
        {
            return static_cast<T>(cur + static_cast<T>(n));
        }

        friend counting_iterator operator+(counting_iterator it, difference_type n)
        {
            it += n;
            return it;
        }

        friend counting_iterator operator+(difference_type n, counting_iterator it)
        {
            it += n;
            return it;
        }

        friend counting_iterator operator-(counting_iterator it, difference_type n)
        {
            it -= n;
            return it;
        }

        friend difference_type operator-(const counting_iterator& a, const counting_iterator& b)
        {
            return static_cast<difference_type>(a.cur) - static_cast<difference_type>(b.cur);
        }

        friend bool operator==(const counting_iterator& a, const counting_iterator& b)
        {
            return a.cur == b.cur;
        }

        friend auto operator<=>(const counting_iterator& a, const counting_iterator& b)
        {
            return a.cur <=> b.cur;
        }
    };
}