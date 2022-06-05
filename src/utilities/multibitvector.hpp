#ifndef MULTIBITVECTOR_HPP
#define MULTIBITVECTOR_HPP
#include "utilities/bitvector.hpp"

template<bool track_count>
class Multibitvector
{
  static constexpr const size_t bitvector_mask = size_t(63);
  static constexpr const int bitvector_log = 6;
  size_t bits,
         set_bits,
         num_bitvectors;
  Bitvector *multibitvector;
  using DivResult = struct
  {
    size_t quot;
    size_t rem;
  };
  DivResult quot_rem_get(size_t idx) const noexcept
  {
    DivResult divresult;
    divresult.quot = idx >> bitvector_log;
    divresult.rem = idx & bitvector_mask;
    return divresult;
  }

  public:
  Multibitvector(size_t _bits)
    : bits(_bits)
    , set_bits(0)
  {
    DivResult divresult = quot_rem_get(bits);
    num_bitvectors = divresult.quot + ((divresult.rem > 0) ? 1 : 0);
    multibitvector = new Bitvector[num_bitvectors];
  }
  Multibitvector(const Multibitvector& m) // Copy constructor
    : bits(m.bits)
    , set_bits(m.set_bits)
    , num_bitvectors(m.num_bitvectors)
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
    DivResult divresult = quot_rem_get(idx);
    if constexpr (track_count)
    {
      set_bits += !(multibitvector[divresult.quot][divresult.rem]);
    }
    multibitvector[divresult.quot].set(divresult.rem);
  }
  void reset(size_t idx) noexcept
  {
    DivResult divresult = quot_rem_get(idx);
    assert(idx < bits &&
           (!multibitvector[divresult.quot][divresult.rem] ||
            set_bits > 0));
    if constexpr (track_count)
    {
      set_bits -= multibitvector[divresult.quot][divresult.rem];
    }
    multibitvector[divresult.quot].reset(divresult.rem);
  }
  bool operator[](size_t idx) const noexcept
  {
    assert(idx < bits);
    DivResult divresult = quot_rem_get(idx);
    return multibitvector[divresult.quot][divresult.rem];
  }
  size_t size(void) const noexcept
  {
    return bits;
  }
  size_t count(void) const noexcept
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
    }
  }
};
#endif
