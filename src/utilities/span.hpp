#ifndef SPAN_HPP
#define SPAN_HPP
#include <cstddef>
#include <cassert>

template<typename T>
class Span
{
   T* ptr;
   size_t num_elements;

public:
    Span(T* _ptr, size_t _num_elements) noexcept
      : ptr(_ptr),
        num_elements(_num_elements)
    {}

    T const& operator[](size_t idx) const noexcept
    {
      assert(idx < num_elements);
      return ptr[idx];
    }

    T const& operator++(void) noexcept
    {
      ptr++;
      return *this;
    }

    [[nodiscard]] size_t size(void) const noexcept { return num_elements; }

    T* begin(void) noexcept
    {
      return ptr;
    }

    T* end(void) noexcept
    {
      return ptr + num_elements;
    }

    T* data(void) noexcept
    {
      return ptr;
    }
};
#endif
