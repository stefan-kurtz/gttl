#ifndef MULTIBITVECTOR_HPP
#define MULTIBITVECTOR_HPP
#include <cstddef>
#include <cassert>
#include "utilities/bitvector.hpp"

template<bool track_count>
class Multibitvector
{
  static constexpr const size_t bitvector_mask = size_t(63);
  static constexpr const int bitvector_log = 6;
  const size_t num_bits,
               num_bitvectors;
  size_t set_bits;
  Bitvector *multibitvector;
  public:
  Multibitvector(size_t _bits)
    : num_bits(_bits)
    , num_bitvectors((num_bits + 63) >> bitvector_log)
    , set_bits(0)
    , multibitvector(new Bitvector [num_bitvectors])
  { }
  Multibitvector(const Multibitvector& m) // Copy constructor
    : num_bits(m.num_bits)
    , num_bitvectors(m.num_bitvectors)
    , set_bits(m.set_bits)
    , multibitvector(new Bitvector[m.num_bitvectors])
  {
    for (size_t idx = 0; idx < num_bitvectors; idx++)
    {
      multibitvector[idx] = m.multibitvector[idx];
    }
  }
  ~Multibitvector(void)
  {
    delete[] multibitvector;
  }
  void set(size_t idx)
  {
    if constexpr (track_count)
    {
      set_bits += not (multibitvector[idx >> bitvector_log]
                                     [idx & bitvector_mask]);
    }
    multibitvector[idx >> bitvector_log].set(idx & bitvector_mask);
  }
  void reset(size_t idx) noexcept
  {
    assert(idx < num_bits);
    if constexpr (track_count)
    {
      set_bits -= multibitvector[idx >> bitvector_log][idx & bitvector_mask];
    }
    multibitvector[idx >> bitvector_log].reset(idx & bitvector_mask);
  }
  bool operator[](size_t idx) const noexcept
  {
    assert(idx < num_bits);
    return multibitvector[idx >> bitvector_log][idx & bitvector_mask];
  }
  [[nodiscard]] size_t size(void) const noexcept { return num_bits; }

  [[nodiscard]] size_t size_in_bytes(void) const
  {
    return num_bitvectors * sizeof(Bitvector);
  }
  void operator |=(const Multibitvector &other) noexcept
  {
    static_assert(not track_count);
    for (size_t idx = 0; idx < num_bitvectors; idx++)
    {
      multibitvector[idx] |= other.multibitvector[idx];
    }
  }
  bool operator==(const Multibitvector<track_count> &rhs) const noexcept
  {
    if (num_bits != rhs.num_bits)
    {
      return false;
    }
    if constexpr (track_count)
    {
      if (set_bits != rhs.set_bits)
      {
        return false;
      }
    }
    for (size_t idx = 0; idx < num_bitvectors; idx++)
    {
      if (multibitvector[idx] != rhs.multibitvector[idx])
      {
        return false;
      }
    }
    return true;
  }
  bool operator!=(const Multibitvector<track_count> &rhs) const noexcept
  {
    return !(*this == rhs);
  }
  [[nodiscard]] size_t count(void) const noexcept
  {
    if constexpr (track_count)
    {
      return set_bits;
    } else
    {
      size_t this_count = 0;
      for (size_t idx = 0; idx < num_bitvectors; idx++)
      {
        this_count += multibitvector[idx].count();
      }
      return this_count;
    }
  }
};
#endif
