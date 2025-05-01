#ifndef FIXED_SIZE_PQ_HPP
#define FIXED_SIZE_PQ_HPP
#include <cassert>
#include <cstdlib>

template <typename BaseType, int (*compare)(const BaseType *, const BaseType *)>
class FixedSizePriorityQueue
{
 private:
  size_t capacity, numofelements;
  BaseType *elements;

 public:
  FixedSizePriorityQueue(size_t maxnumofelements)
  {
    /* one extra element for overflow */
    elements = new BaseType[maxnumofelements + 1];
    capacity = maxnumofelements;
    numofelements = 0;
  }

  ~FixedSizePriorityQueue(void)
  {
    delete[] elements;
  }

  void reset(void) { numofelements = 0; }
  bool is_empty(void) const { return static_cast<bool>(numofelements == 0); }
  bool is_full(void) const
  {
    return static_cast<bool>(numofelements == capacity);
  }
  size_t size(void) const { return numofelements; }

  void add(const BaseType *value)
  {
    BaseType *ptr;

    /* store elements in reverse order, i.e.\ with the element of maximum
       priority at the last index */
    for (ptr = elements + numofelements;
         ptr > elements && compare(value, ptr - 1) < 0; ptr--)
    {
      *ptr = *(ptr - 1);
    }
    if (is_full())
    {
      if (ptr < elements + numofelements)
      {
        *ptr = *value;
      }
    } else
    {
      *ptr = *value;
      numofelements++;
    }
  }

  BaseType &operator[](size_t idx)
  {
    assert(idx < size());
    return elements[idx];
  }
};
#endif
