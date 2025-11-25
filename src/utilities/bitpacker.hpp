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
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <type_traits>
#include "utilities/mathsupport.hpp"

/* class is restricted to the case where a bitgroup_size is <= bits_basetype
   and all, but the last bitgroup fits into the first value of type
   basetype */

template <int sizeof_unit,int bit_groups>
struct GttlBitPacker
{
  private:
    static constexpr const bool show_bit_groups = false;
    std::array<int, bit_groups> bit_group_sizes;
  public:
    std::array<int, bit_groups> shift_tab;
    std::array<uint64_t, bit_groups> mask_tab;
    /* the following two are only needed for
       sizeof_unit > sizeof(basetype) */
    int overflow_left_shift;
    uint64_t max_overflow;

    GttlBitPacker(void) = default;
    GttlBitPacker(const std::array<int, bit_groups> _bit_group_sizes) :
      overflow_left_shift(0),
      max_overflow(0)
    {
      using basetype = std::conditional_t<
                                   sizeof_unit >= 8,
                                   uint64_t,
                                   std::conditional_t<sizeof_unit >= 4,
                                                      uint32_t,
                                                      uint16_t>>;
      static_assert(bit_groups >= 2 &&
                    sizeof_unit >= sizeof(basetype) &&
                    sizeof_unit <= sizeof(basetype) + 7);
      static constexpr const int bits_basetype = CHAR_BIT * sizeof(basetype);
      int count = 0;
      int idx;
      for (idx = 0; idx < bit_groups; idx++)
      {
        bit_group_sizes[idx] = _bit_group_sizes[idx];
        assert(_bit_group_sizes[idx] <= bits_basetype);
        if (count + _bit_group_sizes[idx] > bits_basetype)
        {
          assert(idx < bit_groups);
          break;
        }
        count += _bit_group_sizes[idx];
        shift_tab[idx] = bits_basetype - count;
        if constexpr (show_bit_groups)
        {
          std::cout << "# group " << idx << ", size="  << _bit_group_sizes[idx]
                    << ", shift=" << shift_tab[idx] << '\n';
        }
        mask_tab[idx] = gttl_bits2maxvalue<uint64_t>(_bit_group_sizes[idx]);
      }
      int overflow_bits = 0;
      if (idx < bit_groups)
      {
        assert(idx + 1 == bit_groups && sizeof_unit > sizeof(basetype));
        overflow_bits = _bit_group_sizes[bit_groups-1];
      }
      assert(count <= bits_basetype);
      const int remaining_bits = bits_basetype - count;
      if constexpr (sizeof_unit > sizeof(basetype))
      {
      /* if remaining_bits >= overflow_bits, then overflow value would fit
         into basetype and this does not require an overflow value */
        assert(remaining_bits < overflow_bits);
        assert(static_cast<size_t>(overflow_bits) <= remaining_bits +
                                                     (sizeof_unit -
                                                      sizeof(basetype)) *
                                                      CHAR_BIT);
        assert(overflow_bits <= bits_basetype);
        max_overflow = gttl_bits2maxvalue<uint64_t>(overflow_bits);
        overflow_left_shift = overflow_bits - remaining_bits;
      }
    }

    [[nodiscard]] int bit_group_size_get(int idx) const noexcept
    {
      assert(idx < bit_groups);
      return bit_group_sizes[idx];
    }

    void pretty_print(const char *tag) const noexcept
    {
      std::cout << "# " << tag << '\n';
      for (int idx = 0; idx < bit_groups; idx++)
      {
        std::cout << "# bit_group\t" << idx << "\t" << bit_group_sizes[idx]
                  << '\n';
      }
      std::cout << "# overflow_left_shift\t" << overflow_left_shift
                << '\n';
      std::cout << "# max_overflow\t" << max_overflow << '\n';
    }
};
#endif
