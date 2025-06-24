#ifndef FS_PRIO_STORE_HPP
#define FS_PRIO_STORE_HPP
#include <functional>
#include <vector>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <string>

template <typename BaseType, bool heap_based = true>
class FSPrioStore
{
  std::vector<BaseType> elements;
  const size_t capacity;
  size_t num_elements; /* only used if heap_based is false */
  bool sorted;         /* only used if heap_based is true */
 public:
  FSPrioStore(size_t _capacity)
    : elements({})
    , capacity(_capacity)
    , num_elements(0)
    , sorted(false)
  {
    if constexpr (not heap_based)
    {
      elements.resize(capacity + 1); // need one extra to move unnecessay elem
    }
  }

  void reset(void)
  {
    elements.clear();
    if constexpr (not heap_based)
    {
      num_elements = 0;
    }
  }

  size_t size(void) const
  {
    if constexpr (heap_based)
    {
      return elements.size();
    } else
    {
      return num_elements;
    }
  }

  bool is_empty(void) const
  {
    return static_cast<bool>(size() == 0);
  }

  bool is_full(void) const
  {
    return static_cast<bool>(size() == capacity);
  }

  void add(const BaseType value)
  {
    if constexpr (heap_based)
    {
      if (not is_full())
      {
        elements.push_back(value);
        if (is_full())
        {
          std::make_heap(elements.begin(),elements.end(),
                         std::greater<BaseType>());
        }
      } else
      {
        if (value > elements.front())
        {
          elements.push_back(value);
          std::push_heap(elements.begin(),elements.end(),
                         std::greater<BaseType>());
          assert(std::is_heap(elements.begin(),elements.end() - 1,
                              std::greater<BaseType>()));
          std::pop_heap(elements.begin(),elements.end(),
                        std::greater<BaseType>());
          elements.pop_back();
        }
      }
    } else
    {
      size_t idx;
      for (idx = num_elements; idx > 0 and value > elements[idx-1]; idx--)
      {
        elements[idx] = elements[idx-1];
      }
      if (not is_full() or idx < num_elements)
      {
        elements[idx] = value;
      }
      num_elements += (is_full() ? 0 : 1);
    }
  }
  void sort(void)
  {
    if constexpr (heap_based)
    {
      if (not sorted)
      {
        if (not is_heap(elements.begin(),elements.end(),
                        std::greater<BaseType>()))
        {
          std::make_heap(elements.begin(),elements.end(),
                         std::greater<BaseType>());
        }
        std::sort_heap(elements.begin(),elements.end(),
                       std::greater<BaseType>());
        sorted = true;
      }
    }
  }
  using Iterator = typename std::vector<BaseType>::iterator;
  using ConstIterator = typename std::vector<BaseType>::const_iterator;
  Iterator begin()
  {
    if constexpr (heap_based)
    {
      assert(sorted);
    }
    return elements.begin();
  }
  Iterator end()
  {
    if constexpr (heap_based)
    {
      return elements.end();
    } else
    {
      return elements.begin() + num_elements;
    }
  }
  ConstIterator begin() const
  {
    if constexpr (heap_based)
    {
      assert(sorted);
    }
    return elements.begin();
  }
  ConstIterator end() const
  {
    if constexpr (heap_based)
    {
      return elements.end();
    } else
    {
      return elements.begin() + num_elements;
    }
  }
  ConstIterator cbegin() const
  {
    if constexpr (heap_based)
    {
      assert(sorted);
    }
    return elements.cbegin();
  }
  ConstIterator cend() const
  {
    if constexpr (heap_based)
    {
      return elements.cend();
    } else
    {
      return elements.cbegin() + num_elements;
    }
  }
  void show(void) const
  {
    if constexpr (heap_based)
    {
      assert(sorted);
    }
    std::cout << "[";
    for (size_t j = 0; j < size(); j++)
    {
      if (j > 0)
      {
        std::cout << ",";
      }
      std::cout << std::to_string(elements[j]);
    }
    std::cout << "]\n";
  }
};
#endif
