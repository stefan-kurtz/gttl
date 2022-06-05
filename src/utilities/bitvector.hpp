#ifndef BITVECTOR_HPP
#define BITVECTOR_HPP
#include <cassert>
#include <climits>
#include <cstdint>
#include <string>
#include "utilities/mathsupport.hpp"

class Bitvector
{
  static constexpr const size_t bits = sizeof(uint64_t) * CHAR_BIT;
  static constexpr const uint64_t first_bit = uint64_t(1) << (bits - 1);
  uint64_t value;
  uint64_t the_bit(size_t idx) const noexcept
  {
    assert(idx < bits);
    return first_bit >> idx;
  }
  public:
  Bitvector(void)
    : value(0)
    {}
  Bitvector(uint64_t _value)
    : value(_value)
    {}
  void set(size_t idx)
  {
    value |= the_bit(idx);
  }
  void reset(size_t idx)
  {
    value &= ~the_bit(idx);
  }
  bool operator[](size_t idx) const noexcept
  {
    return static_cast<bool>(value & the_bit(idx));
  }
  size_t count(void) const noexcept
  {
    return popcount_uint64_t(value);
  }
  std::string to_string(void) const noexcept
  {
    std::string s{};
    s.reserve(bits);
    for (size_t idx = 0; idx < bits; idx++)
    {
      s += (*this)[idx] ? '1' : '0';
    }
    return s;
  }
};
#endif
