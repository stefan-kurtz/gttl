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

    T const& operator++() noexcept
    {
      ptr++;
      return *this;
    }

    size_t size() const noexcept
    {
      return num_elements;
    }

    T* begin() noexcept
    {
      return ptr;
    }

    T* end() noexcept
    {
      return ptr + num_elements;
    }
};
#endif
