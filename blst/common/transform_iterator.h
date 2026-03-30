#include <iterator>
#include <type_traits>

namespace blst
{

    template<class Iter, class Unary>
    struct transform_iterator
    {
        using value_type = std::invoke_result_t<Unary, typename std::iterator_traits<Iter>::value_type>;
        using difference_type = typename std::iterator_traits<Iter>::difference_type;
        using iterator_category = std::random_access_iterator_tag;
        using reference = value_type;
        using pointer = void;

        Iter it;
        Unary fun;

        transform_iterator() = default;
        BLST_HD transform_iterator(Iter i, Unary f) : it(i), fun(f) {}

        BLST_HD reference operator*() const
        {
            return fun(*it);
        }

        BLST_HD transform_iterator& operator++()
        {
            ++it;
            return *this;
        }

        BLST_HD transform_iterator operator++(int)
        {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }

        BLST_HD transform_iterator& operator--()
        {
            --it;
            return *this;
        }

        BLST_HD transform_iterator operator--(int)
        {
            auto tmp = *this;
            --(*this);
            return tmp;
        }

        BLST_HD transform_iterator& operator+=(difference_type n)
        {
            it += n;
            return *this;
        }

        BLST_HD transform_iterator& operator-=(difference_type n)
        {
            it -= n;
            return *this;
        }

        BLST_HD reference operator[](difference_type n) const
        {
            return fun(it[n]);
        }

        BLST_HD friend transform_iterator operator+(transform_iterator t, difference_type n)
        {
            t += n;
            return t;
        }

        BLST_HD friend transform_iterator operator+(difference_type n, transform_iterator t)
        {
            t += n;
            return t;
        }

        BLST_HD friend transform_iterator operator-(transform_iterator t, difference_type n)
        {
            t -= n;
            return t;
        }

        BLST_HD friend difference_type operator-(const transform_iterator& a, const transform_iterator& b)
        {
            return a.it - b.it;
        }

        BLST_HD friend bool operator==(const transform_iterator& a, const transform_iterator& b)
        {
            return a.it == b.it;
        }

        BLST_HD friend auto operator<=>(const transform_iterator& a, const transform_iterator& b)
        {
            return a.it <=> b.it;
        }
    };

}