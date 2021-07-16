#ifndef FS_PRIO_STORE_HPP
#define FS_PRIO_STORE_HPP
#include <cassert>
#include <cstdbool>
#include <cstdlib>

template<typename T>
using CompareFunc = bool (*)(const T *, const T *);

template <typename T>
class FSPrioStore
{
 private:
  T *elements;
  size_t capacity, numofelements;

 public:
  FSPrioStore(size_t _capacity)
  {
    /* one extra element for overflow */
    elements = new T[_capacity + 1];
    capacity = _capacity;
    numofelements = 0;
  }

  ~FSPrioStore(void)
  {
    delete[] elements;
  }

  void reset(void)
  {
    numofelements = 0;
  }

  bool is_empty(void) const
  {
    return static_cast<bool>(numofelements == 0);
  }

  bool is_full(void) const
  {
    return static_cast<bool>(numofelements == capacity);
  }

  size_t size(void) const
  {
    return numofelements;
  }

  void add(const T *value,CompareFunc<T> compare)
  {
    T *ptr;

    for (ptr = elements + numofelements;
         ptr > elements && compare(value, ptr - 1); ptr--)
    {
      *ptr = *(ptr - 1);
    }
    if (!is_full() || ptr < elements + numofelements)
    {
      *ptr = *value;
    }
    numofelements += (this->is_full() ? 0 : 1);
  }

  T &operator[](size_t idx) const
  {
    assert(idx < size());
    return elements[idx];
  }
};
#endif
