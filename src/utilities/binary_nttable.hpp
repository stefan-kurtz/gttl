#ifndef BINARY_NTTABLE_HPP
#define BINARY_NTTABLE_HPP

#include <cmath>
#include <cstddef>
#include <cstdint>
#include "utilities/multibitvector.hpp"

class BinaryNtTable
{
  private:
  const size_t s_value;
  const size_t r_value;
  const uint64_t s_mask;
  const uint64_t r_mask;
  Multibitvector<false> table;

  public:
  BinaryNtTable(size_t _s, size_t _r)
    : s_value(_s)
    , r_value(_r)
    , s_mask(_s == 0 ? 0 : ((~uint64_t(0)) << (64 - _s)))
    , r_mask(               (~uint64_t(0)) >> (64 - _r))
    , table(Multibitvector<false>(size_t(1) << _r))
  { }

  void merge(const BinaryNtTable &other)
  {
    assert(s_value == other.s_value and r_value == other.r_value);
    table |= other.table;
  }

  void add_hash(uint64_t hash)
  {
    if ((hash & s_mask) != 0)
    {
      return;
    }

    table.set(static_cast<size_t>(hash & r_mask));
  }

  double estimate_F0(void) const
  {
    const size_t p0 = (size_t(1) << r_value) - table.count();

    if (p0 == 0)
    {
      throw std::string("a count of zero leads to undefined values");
    }
    return -std::log(static_cast<double>(p0) / (uint64_t(1) << r_value)) *
           (uint64_t(1) << (s_value + r_value));
  }
};
#endif