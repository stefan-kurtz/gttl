/*
  Copyright (c) 2021 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2021 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/
#ifndef MATHSUPPORT_HPP
#define MATHSUPPORT_HPP
#include <cstdint>
#include <cstddef>
#include <climits>
#include <cassert>
#include <cmath>
#include <numeric>

#ifndef __has_builtin         // Optional of course.
#define __has_builtin(X) 0  // Compatibility with non-clang compilers.
#endif

template<typename T,int bits>
inline constexpr T gttl_bits2maxvalue(void)
{
  static_assert(static_cast<T>(-1) >= 0);
  if constexpr (bits == CHAR_BIT * sizeof(T))
  {
    return ~static_cast<T>(0);
  }
  static_assert(static_cast<size_t>(bits) < CHAR_BIT * sizeof(T));
  return (static_cast<T>(1) << bits) - 1;
}

template<typename T>
inline T gttl_bits2maxvalue(int bits)
{
  static_assert(static_cast<T>(-1) >= 0);
  if (bits == CHAR_BIT * sizeof(T))
  {
    return ~static_cast<T>(0);
  }
  assert(static_cast<size_t>(bits) < CHAR_BIT * sizeof(T));
  return (static_cast<T>(1) << bits) - 1;
}

template<typename T>
inline T gttl_bits2maxvalue_not_full(int bits)
{
  static_assert(static_cast<T>(-1) >= 0);
  assert(static_cast<size_t>(bits) < CHAR_BIT * sizeof(T));
  return (static_cast<T>(1) << bits) - 1;
}

template<typename Numtype>
inline int gttl_required_bits(Numtype value)
{
  if (value == 0)
  {
    return 0;
  }
#if __has_builtin(__builtin_clz)
#include <climits>
  static_assert(sizeof(value) <= sizeof(unsigned long));
  return sizeof(unsigned long) * CHAR_BIT -
         __builtin_clzl(static_cast<unsigned long>(value));
#else
  int count;
  for(count = 0; value > 0; count++)
  {
    value >>= 1;
  }
  return count;
#endif
}

inline size_t popcount_uint64_t(uint64_t value)
{
#if __has_builtin(__builtin_popcountl)
  return __builtin_popcountl(static_cast<long>(value));
#else
  size_t pc = 0;
  for (; value != 0; value &= value - 1)
  {
    pc++;
  }
  return pc;
#endif
}

inline double mega_bytes(size_t bytes)
{
  return static_cast<double>(bytes)/(size_t(1024) * size_t(1024));
}

/* compute base 2 logarithm at compile time:
   use as Log2_CT<64>::VALUE; */

template <int numerus>
struct Log2_CT
{
  enum { VALUE = Log2_CT<(numerus+1)/2>::VALUE + 1 };
};
template < > struct Log2_CT<1> { enum { VALUE = 0 }; };
template < > struct Log2_CT<0> { enum { VALUE = 0 }; };


/* compute power of base 2 at compile time
   use as Pow_CT<16>::VALUE; */

template <size_t base,int numerus>
struct Pow_CT
{
  enum { VALUE = Pow_CT<base,numerus-1>::VALUE * base };
};
template <size_t base> struct Pow_CT<base,1> { enum { VALUE = base }; };
template <size_t base> struct Pow_CT<base,0> { enum { VALUE = 1 }; };

inline double error_percentage_get(size_t distance,size_t aligned_len)
{
  if (aligned_len == 0)
  {
    assert(distance == 0);
    return 0;
  }
  return 100.0 * static_cast<double>(distance)/(aligned_len/2.0);
}

template<class InputIt>
double gttl_variance(InputIt first,InputIt last)
{
  const auto sum = std::accumulate(first,last,0);
  const auto num_elements = last - first;
  const double mean = static_cast<double>(sum)/num_elements;
  double squared_difference = 0;
  for (auto it = first; it != last; ++it)
  {
    const double diff = *it - mean;
    squared_difference += (diff * diff);
  }
  return squared_difference / num_elements;
}

template<class InputIt>
double gttl_stddev(InputIt first,InputIt last)
{
  return std::sqrt(gttl_variance(first,last));
}
#endif
