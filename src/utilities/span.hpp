#include <cstddef>

template<typename T>
class Span {
   T* ptr;
   size_t len;

public:
    Span(T* _ptr, size_t _len) noexcept
        : ptr(_ptr),
          len(_len)
    {}

    T const& operator[](size_t idx) const noexcept {
        return ptr[idx];
    }

    size_t size() const noexcept {
        return len;
    }

    T* begin() noexcept {
        return ptr;
    }

    T* end() noexcept {
        return ptr + len;
    }
};
