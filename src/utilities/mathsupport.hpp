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

#ifndef __has_builtin         // Optional of course.
#define __has_builtin(X) 0  // Compatibility with non-clang compilers.
#endif

template<typename Numtype>
inline int gt_required_bits(Numtype value)
{
  if (value == 0)
  {
    return 0;
  }
#if __has_builtin(__builtin_clz)
#include <climits>
  static_assert(sizeof(value) <= sizeof(unsigned long));
  return sizeof(value) * CHAR_BIT - __builtin_clzl((unsigned long) value);
#else
  int count;
  for(count = 0; value > 0; count++)
  {
    value >>= 1;
  }
  return count;
#endif
}

inline double mega_bytes(size_t bytes)
{
  return (double) bytes/(static_cast<size_t>(1024 * 1024));
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
#endif
