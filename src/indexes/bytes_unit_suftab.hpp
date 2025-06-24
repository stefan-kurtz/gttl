/*
  Copyright (c) 2022 Stefan Kurtz <stefan.kurtz@uni-hamburg.de>
  Copyright (c) 2022 Center for Bioinformatics, University of Hamburg

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

#ifndef BYTES_UNIT_SUFTAB_HPP
#define BYTES_UNIT_SUFTAB_HPP
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>
#include <climits>
#include "utilities/bitpacker.hpp"
#include "utilities/bytes_unit.hpp"
#include "indexes/gttl_suffixarray.hpp"

template <int sizeof_unit>
class BytesUnitSUFtab
{
  private:
  const size_t nonspecial_suffixes;
  const int sequences_number_bits,
            sequences_length_bits,
            first_group_bits;
  const GttlBitPacker<sizeof_unit,2> bp;
  const BytesUnit<sizeof_unit,2> *bu_suftab;

  std::pair<uint32_t,uint32_t> decode_seqnum_relpos(
                                  const BytesUnit<sizeof_unit,2> &suffix)
                                    const noexcept
  {
    const uint64_t leafnumber_seqnum = suffix.template decode_at<0>(bp);
    const uint64_t leafnumber_relpos = suffix.template decode_at<1>(bp);
    assert(leafnumber_relpos <= UINT32_MAX && leafnumber_seqnum <= UINT32_MAX);
    return {static_cast<uint32_t>(leafnumber_seqnum),
            static_cast<uint32_t>(leafnumber_relpos)};
  }
  const BytesUnit<sizeof_unit,2> *bu_suftab_get(bool with_mmap,
                                                const GttlSuffixArray
                                                  *suffixarray)
    const
  {
    if (with_mmap)
    {
      const uint8_t* suftab_bytes = suffixarray->get_mmap_suftab_bytes();
      return reinterpret_cast<const BytesUnit<sizeof_unit,2> *>(suftab_bytes);
    }
    const std::vector<uint8_t> &suftab_bytes = suffixarray->get_suftab_bytes();
    return reinterpret_cast<const BytesUnit<sizeof_unit,2> *>
                           (suftab_bytes.data());
  }
  public:
  BytesUnitSUFtab(bool with_mmap,
                  const GttlSuffixArray *suffixarray,
                  size_t _nonspecial_suffixes)
    : nonspecial_suffixes(_nonspecial_suffixes)
    , sequences_number_bits(suffixarray->sequences_number_bits_get())
    , sequences_length_bits(suffixarray->sequences_length_bits_get())
    , first_group_bits(suffixarray->sequences_number_get() == 1
                         ? sizeof_unit * CHAR_BIT - sequences_length_bits
                         : sequences_number_bits)
    , bp({first_group_bits,sequences_length_bits})
    , bu_suftab(bu_suftab_get(with_mmap,suffixarray))
  {
    assert(sequences_number_bits + sequences_length_bits <=
           sizeof_unit * CHAR_BIT);
    assert(suffixarray->sequences_number_get() > 1 ||
           sequences_number_bits == 0);
    assert(suffixarray->sequences_number_get() == 1 ||
           sequences_number_bits > 0);
  }
  std::pair<uint32_t,uint32_t> operator[](size_t interval_bound) const noexcept
  {
    return decode_seqnum_relpos(bu_suftab[interval_bound]);
  }
  size_t size(void) const noexcept
  {
    return nonspecial_suffixes;
  }
};
#endif
