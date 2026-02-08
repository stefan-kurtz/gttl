#ifndef MYERSAPM_HPP
#define MYERSAPM_HPP
#include <cstddef>
#include <cstdlib>
#include <cassert>
#include <climits>
#include <format>
#include <array>
#include <stdexcept>

template<typename T>
class MyersBitvectorAlgorithm
{
  const size_t pattern_length;
  const T bitmask;
  T Mv,
    Pv; /* 1^{m} */
  std::array<T,UCHAR_MAX+1> eq_bits;
  size_t cost;
  public:
  void reset_for_next_sequence(void)
  {
    Mv = static_cast<T>(0);
    Pv = ~static_cast<T>(0); /* 1^{m} */
    cost = pattern_length;
  }
  MyersBitvectorAlgorithm(const std::string &pattern)
    : pattern_length(pattern.size())
    , bitmask(static_cast<T>(1) << (pattern_length - 1))  /* 10^{m-1} */
  {
    constexpr const size_t max_pattern_length = sizeof(T) * CHAR_BIT;
    if (pattern_length > max_pattern_length)
    {
      throw std::format("MyersBitvectorAlgorithm cannot handle patterns of "
                        "length {}; the maximum pattern length is {}" ,
                        pattern_length, max_pattern_length);
    }
    eq_bits.fill(0);
    for (size_t idx = 0; idx < pattern_length; idx++)
    {
      eq_bits[(int) pattern[idx]] |= (static_cast<T>(1) << idx);
    }
    reset_for_next_sequence();
  }
  size_t transform(char cc)
  {
    const T Eq = eq_bits[(int) cc];
    const T Xh = Eq | (((Eq & Pv) + Pv) ^ Pv);
    T Ph = Mv | ~ (Xh | Pv);
    const T Mh = Pv & Xh;
    cost = cost + ((Ph & bitmask) ? size_t(1) : 0)
                - ((Mh & bitmask) ? size_t(1) : 0);
    const T Xv = Eq | Mv;
    Ph <<= 1; /* shift only once and use value twice */
    Pv = (Mh << 1) | ~ (Xv | Ph);
    Mv = Ph & Xv;
    return cost;
  }
};
#endif
