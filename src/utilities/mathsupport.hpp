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
#include <cstddef>
#include <climits>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

#ifndef __has_builtin         // Optional of course.
#define __has_builtin(X) 0  // Compatibility with non-clang compilers.
#endif

template<typename T,int bits>
inline consteval T gttl_bits2maxvalue(void)
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

inline constexpr const double bytesPerMegaByte = size_t{1024} * size_t{1024};

inline constexpr double mega_bytes(size_t bytes) noexcept
{
  return static_cast<double>(bytes)/bytesPerMegaByte;
}

template<typename T>
static constexpr T gttl_safe_power(T a, T b)
{
  static_assert(std::is_integral_v<T> and std::is_unsigned_v<T>);

  T prod = 1;
  constexpr const T max_value = std::numeric_limits<T>::max();

  while (b > 0)
  {
    if ((b & 1) == 1)
    {
      if (a != 0 and prod > max_value / a)
      {
        throw std::overflow_error(
              std::string("overflow when evaluating ")
            + std::to_string(prod) + " * " + std::to_string(a)
            + std::string(" to compute pow(")
            + std::to_string(a) + std::string(", ")
            + std::to_string(b) + std::string(")"));
      }
      prod *= a;
    }
    b >>= 1;
    if (b > 0 and a != 0 and a > max_value / a)
    {
        throw std::overflow_error(
              std::string("overflow when evaluating ")
            + std::to_string(prod) + " * " + std::to_string(a)
            + std::string(" to compute pow(")
            + std::to_string(a) + std::string(", ")
            + std::to_string(b) + std::string(")"));
    }
    a *= a;
  }
  return prod;
}

inline double error_percentage_get(size_t distance,size_t aligned_len)
{
  if (aligned_len == 0)
  {
    assert(distance == 0);
    return 0;
  }
  return 100.0 * static_cast<double>(distance)/(aligned_len/2.0);
}

// Welford's Algorithm for numerically determining the variance of a population
template<class InputIt>
double gttl_variance(InputIt first, InputIt last)
{
  size_t n = 0;
  double mean = 0.0;
  double M2 = 0.0;
  for(; first != last; ++first)
  {
    ++n;
    const double x = static_cast<double>(*first);
    const double delta = x - mean;
    mean += delta / n;
    M2 += delta * (x - mean);
  }
  return (n > 0) ? M2 / n : 0.0;
}

template<class InputIt>
double gttl_stddev(InputIt first,InputIt last)
{
  return std::sqrt(gttl_variance(first,last));
}
#endif
