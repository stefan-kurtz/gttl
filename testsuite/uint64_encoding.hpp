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
#ifndef UINT64_ENCODING_HPP
#define UINT64_ENCODING_HPP

#include <cassert>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>

/* Class Uint64Encoding
 * handling the encoding and decoding of given numbers in one uint64_t (64bits)
 */

template <int bit_groups>
class Uint64Encoding
{
  private:
    std::array<int, bit_groups> bit_group_sizes;
    std::array<int, bit_groups> shift_tab;
    std::array<uint64_t, bit_groups> mask_tab;

  public:
   /* Constructor takes array, the number of values in array will be the count
      of numbers in code, the values are for the number of bits which will
      be used for the values in en/decoding, make shure it summs to not more
      than 64 bits */
    Uint64Encoding(void) {}
    Uint64Encoding(const std::array<int, bit_groups> _bit_group_sizes)
    {
      int count = 0;
      for (int idx = 0; idx < bit_groups; idx++)
      {
        bit_group_sizes[idx] = _bit_group_sizes[idx];
        assert(_bit_group_sizes[idx] < 64);
        count += _bit_group_sizes[idx];
        assert(count <= 64);
        shift_tab[idx] = 64 - count;
        mask_tab[idx] = (uint64_t(1) << _bit_group_sizes[idx]) - 1;
      }
    }

    /* Returns the encoded version of the numbers in the given array
    * checks if numbers fitting given bitranges
    */

    [[nodiscard]] uint64_t
    encode(const std::array<uint64_t, bit_groups> &arr) const noexcept
    {
      uint64_t code = 0;

      for (int idx = 0; idx < bit_groups; idx++)
      {
        assert(arr[idx] <= mask_tab[idx]);
        code |= (arr[idx] << shift_tab[idx]);
      }
      return code;
    }

    template <int idx>
    [[nodiscard]] uint64_t decode_at(uint64_t code) const noexcept
    {
      static_assert(idx < bit_groups);
      return (code >> shift_tab[idx]) & mask_tab[idx];
    }

    template <int idx>
    [[nodiscard]] int bit_group_size(void) const noexcept
    {
      static_assert(idx < bit_groups);
      return bit_group_sizes[idx];
    }

    void pretty_print(const char *tag) const noexcept
    {
      std::cout << "# " << tag << '\n';
      for (size_t idx = 0; idx < bit_groups; idx++)
      {
        std::cout << "# bit_group\t" << idx << "\t" << bit_group_sizes[idx]
                  << std::endl;
      }
    }
};
#endif
