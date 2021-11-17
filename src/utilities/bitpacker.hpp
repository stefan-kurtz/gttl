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
#ifndef BITPACKER_HPP
#define BITPACKER_HPP

#include <cassert>
#include <climits>
#include <array>
#include <cstdint>
#include <iostream>
#include "utilities/mathsupport.hpp"

template <typename basetype,int sizeof_unit,int bit_groups>
struct GttlBitPacker
{
  private:
    static_assert(sizeof_unit == sizeof(basetype) ||
                  sizeof_unit == sizeof(basetype) + 1);
    static constexpr const bool show_bit_groups = false;
    std::array<int, bit_groups> bit_group_sizes;
    static constexpr const int bits_basetype = CHAR_BIT * sizeof(basetype);
  public:
    std::array<int, bit_groups> shift_tab;
    std::array<basetype, bit_groups> mask_tab;
    /* the following two are only needed for
       sizeof_unit = sizeof(basetype) + 1 */
    int overflow_left_shift;
    uint64_t max_overflow;
    GttlBitPacker(void) {}
    GttlBitPacker(const std::array<int, bit_groups> _bit_group_sizes) :
      overflow_left_shift(0),
      max_overflow(0)
    {
      int count = 0, idx;
      static_assert(bit_groups >= 2);
      for (idx = 0; idx < bit_groups; idx++)
      {
        if constexpr (show_bit_groups)
        {
          std::cout << "group " << idx << ", size="  << _bit_group_sizes[idx]
                    << std::endl;
        }
        bit_group_sizes[idx] = _bit_group_sizes[idx];
        assert(_bit_group_sizes[idx] <= bits_basetype);
        if (count + _bit_group_sizes[idx] > bits_basetype)
        {
          break;
        }
        count += _bit_group_sizes[idx];
        shift_tab[idx] = bits_basetype - count;
        mask_tab[idx] = gttl_bits2maxvalue<basetype>(_bit_group_sizes[idx]);
      }
      int overflow_bits = 0;
      if (idx < bit_groups)
      {
        assert(sizeof_unit == sizeof(basetype)+1);
        assert(idx+1 == bit_groups);
        overflow_bits = _bit_group_sizes[idx];
      }
      assert(count <= bits_basetype);
      const int remaining_bits = bits_basetype - count;
      if constexpr (sizeof_unit == sizeof(basetype)+1)
      {
      /* if remaining_bits >= overflow_bits, then overflow value would fit
         into basetype and this does not require an overflow value */
        assert(remaining_bits < overflow_bits);
        assert(overflow_bits <= remaining_bits + CHAR_BIT);
        assert(overflow_bits <= bits_basetype);
        max_overflow = gttl_bits2maxvalue<uint64_t>(overflow_bits);
        overflow_left_shift = overflow_bits - remaining_bits;
      }
    }

    int bit_group_size_get(int idx) const noexcept
    {
      assert(idx < bit_groups);
      return bit_group_sizes[idx];
    }

    void pretty_print(const char *tag) const noexcept
    {
      std::cout << "# " << tag << std::endl;
      for (size_t idx = 0; idx < bit_groups; idx++)
      {
        std::cout << "# bit_group\t" << idx << "\t" << bit_group_sizes[idx]
                  << std::endl;
      }
    }
};
#endif
