#ifndef BITVECTOR_HPP
#define BITVECTOR_HPP
#include <cassert>
#include <climits>
#include <cstdint>

class Bitvector
{
  static constexpr const size_t bits = sizeof(uint64_t) * CHAR_BIT;
  uint64_t value;
  uint64_t the_bit(size_t idx) const noexcept
  {
    assert(idx < bits);
    return uint64_t(1) << (bits - 1 - idx);
  }
  public:
  Bitvector(void) :
    value(0) {}
  Bitvector(uint64_t _value) :
    value(_value) {}
  void set(size_t idx)
  {
    value |= the_bit(idx);
  }
  void unset(size_t idx)
  {
    value &= ~the_bit(idx);
  }
  bool is_set(size_t idx) const noexcept
  {
    return static_cast<bool>(value & the_bit(idx));
  }
};
#endif
