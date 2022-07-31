#ifndef SORTED_MATCH_LIST_HPP
#define SORTED_MATCH_LIST_HPP
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
#include <cstddef>
#include <cassert>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <string>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <array>
#include <cmath>
#include <type_traits>
#include <iomanip>
#include <ios>
#include <tuple>

#include "utilities/mathsupport.hpp"
#include "utilities/is_big_endian.hpp"
#include "utilities/unused.hpp"
#include "utilities/ska_lsb_radix_sort.hpp"
#include "utilities/remove_duplicates.hpp"
#include "utilities/bytes_unit.hpp"
#include "sequences/gttl_multiseq.hpp"
#include "sequences/matching_characters.hpp"

template<bool check_bounds,bool (*matching_method)(char,char)>
static inline std::pair<size_t,size_t> maximize_on_both_ends(
                                               const char *seq0,
                                               size_t len0,
                                               size_t abs_start_pos0,
                                               const char *seq1,
                                               size_t len1,
                                               size_t abs_start_pos1,
                                               size_t seed_size)
{
  static auto pointers_leq = [] (const char *a,const char *b)
  {
    if constexpr (check_bounds)
    {
      return a <= b;
    }
    return true;
  };

  assert(abs_start_pos0 + seed_size - 1 < len0 &&
         abs_start_pos1 + seed_size - 1 < len1);
  const char *bwd0, *bwd1;
  for (bwd0 = seq0 + abs_start_pos0 - 1,
       bwd1 = seq1 + abs_start_pos1 - 1;
       pointers_leq(seq0,bwd0) &&
       pointers_leq(seq1,bwd1) &&
       matching_method(*bwd0,*bwd1);
       bwd0--, bwd1--)
       /* Nothing */ ;
  assert(seq0 + abs_start_pos0 >= (bwd0+1));
  const size_t left_extend
    = static_cast<size_t>(seq0 + abs_start_pos0 - 1 - (bwd0+1) + 1);

  const char *fwd0, *fwd1;
  for (fwd0 = seq0 + abs_start_pos0 + seed_size,
       fwd1 = seq1 + abs_start_pos1 + seed_size;
       pointers_leq(fwd0,seq0 + len0 - 1) &&
       pointers_leq(fwd1,seq1 + len1 - 1) &&
       matching_method(*fwd0,*fwd1);
       fwd0++, fwd1++)
         /* Nothing */;
  assert(fwd0 >= seq0 + abs_start_pos0 + seed_size);
  const size_t right_extend
    = static_cast<size_t>(fwd0 - (seq0 + abs_start_pos0 + seed_size));
  return std::make_pair(left_extend,right_extend);
}

#define MAXIMIZE_ON_BOTH_ENDS(CHECK_BOUNDS,MATCHING_CHARACTERS)\
        std::tie(left_extend,right_extend) \
          = maximize_on_both_ends<CHECK_BOUNDS,MATCHING_CHARACTERS>\
                                 (seq0,\
                                  len0,\
                                  pp.startpos0,\
                                  seq1,\
                                  len1,\
                                  pp.startpos1,\
                                  qgram_length)

template<class SeedEnumeratorClass,
         bool self_match,
         int sizeof_unit_match,
         bool seed_output,
         int ref_idx>
class SortedMatchList
{
  static_assert(sizeof_unit_match == 8 || sizeof_unit_match == 9);
  static_assert(ref_idx == 0 || ref_idx == 1);
  static constexpr const int query_idx = ref_idx == 0 ? 1 : 0;
  static constexpr const int ref_pos_idx = ref_idx == 0 ? 2 : 3;
  static constexpr const int query_pos_idx = ref_idx == 0 ? 3 : 2;
  using BytesUnitMatch = BytesUnit<sizeof_unit_match,5>;
  private:
  std::vector<BytesUnitMatch> encoded_match_list;
  size_t minimum_mem_length,
         number_of_seeds,
         number_of_all_matches;
  int bits_for_sequences;
  int remaining_bits_for_length;
  uint64_t maximum_storable_match_length;
  const GttlMultiseq *ref_multiseq,
                     *query_multiseq;
  GttlBitPacker<sizeof_unit_match,5> match_packer;

  public:
  SortedMatchList(size_t qgram_length,
                  size_t _minimum_mem_length,
                  const SeedEnumeratorClass &seed_enumerator,
                  const GttlMultiseq *_ref_multiseq,
                  const GttlMultiseq *_query_multiseq)
    : encoded_match_list({})
    , minimum_mem_length(_minimum_mem_length)
    , number_of_seeds(0)
    , number_of_all_matches(0)
                         /* this must deliver sum of number of bits needed for
                            the sequence positions and the sequence numbers.
                            If all positions refer to the same sequences,
                            then the sequence number are represented by 0
                            bits. */
    , bits_for_sequences(seed_enumerator.sequences_bits_sum_get())
    , remaining_bits_for_length(sizeof_unit_match * CHAR_BIT -
                                bits_for_sequences)
    , maximum_storable_match_length(gttl_bits2maxvalue<uint64_t>
                                       (remaining_bits_for_length))
    , ref_multiseq(_ref_multiseq)
    , query_multiseq(_query_multiseq)
    , match_packer(/* Deliver a std::array<int,5> such that at
                      pos 0: number of bits for the reference sequence numbers
                      pos 1: number of bits for the query_sequence numbers
                      pos 2: number of bits for the reference sequence positions
                      pos 3: number of bits for the query sequence positions
                      pos 4: remaining_bits_for_length.
                      if ref_idx is not 0, then pos 2 and 3 are swapped */
                   seed_enumerator
                   .template match_packer_order_units<ref_idx>
                                                     (remaining_bits_for_length)
                  )
  {
    assert(minimum_mem_length >= qgram_length);
    const size_t length_threshold = minimum_mem_length - qgram_length;

    /* For the following loop to work the
       seed_enumerator must provide iterators begin() and end() which
       when dereferenced deliver an instance of
       struct SortedMatchListPositionPair, or some struct which is compatible.
    */
    for (auto const &pp : seed_enumerator)
    {
      if constexpr (seed_output)
      {
        std::cout << "# seed\t" /* only for debugging and thus not locked */
                  << pp.seqnum0 << "\t"
                  << pp.startpos0 << "\t"
                  << pp.seqnum1 << "\t"
                  << pp.startpos1 << std::endl;
      }
      const char *seq0 = ref_multiseq->sequence_ptr_get(pp.seqnum0);
      const char *seq1 = query_multiseq->sequence_ptr_get(pp.seqnum1);
      const size_t len0 = ref_multiseq->sequence_length_get(pp.seqnum0);
      const size_t len1 = query_multiseq->sequence_length_get(pp.seqnum1);
      size_t left_extend, right_extend;

      if constexpr (self_match)
      {
        if (pp.seqnum0 != pp.seqnum1)
        {
          constexpr const bool check_bounds = true;
          MAXIMIZE_ON_BOTH_ENDS(check_bounds,matching_characters_wc);
        } else
        {
          constexpr const bool check_bounds = false;
          MAXIMIZE_ON_BOTH_ENDS(check_bounds,matching_characters_wc);
        }
      } else
      {
        constexpr const bool check_bounds = false;
        MAXIMIZE_ON_BOTH_ENDS(check_bounds,matching_characters);
      }
      bool seed_okay = false;
      if (left_extend + right_extend >= length_threshold)
      {
        if constexpr (self_match)
        {
          seed_okay = (count_mismatches<matching_characters_wc>
                                       (seq0 + pp.startpos0,
                                        seq1 + pp.startpos1,
                                        qgram_length) == 0);
        } else
        {
          seed_okay = (count_mismatches<matching_characters>
                                       (seq0 + pp.startpos0,
                                        seq1 + pp.startpos1,
                                        qgram_length) == 0);
        }
      }
      if (seed_okay)
      {
        assert(pp.startpos0 >= left_extend && pp.startpos1 >= left_extend);
        const size_t this_match_length = left_extend + qgram_length +
                                         right_extend;
        assert(this_match_length >= minimum_mem_length);
        uint64_t length_stored;
        if (this_match_length >
            minimum_mem_length + maximum_storable_match_length)
        {
          length_stored = maximum_storable_match_length;
          StrFormat msg("cannot store match of length %lu, resort to using "
                        "more space for matches",this_match_length);
          throw msg.str();
        } else
        {
          length_stored = this_match_length - minimum_mem_length;
        }
        if constexpr (ref_idx == 0)
        {
          BytesUnitMatch
            encoded_match(match_packer, /* create match */
                          {pp.seqnum0,
                           pp.seqnum1,
                           pp.startpos0 - left_extend + this_match_length - 1,
                           pp.startpos1 - left_extend + this_match_length - 1,
                           length_stored});
          encoded_match_list.emplace_back(encoded_match);
        } else
        {
          static_assert(ref_idx == 1);
          BytesUnitMatch
            encoded_match(match_packer, /* create match */
                          {pp.seqnum1,
                           pp.seqnum0,
                           pp.startpos1 - left_extend + this_match_length - 1,
                           pp.startpos0 - left_extend + this_match_length - 1,
                           length_stored});
          encoded_match_list.emplace_back(encoded_match);
        }
      }
      number_of_seeds++;
    }
    number_of_all_matches = encoded_match_list.size();
    const bool reversed_byte_order = is_big_endian() ? false : true;
    ska_large_lsb_small_radix_sort(sizeof_unit_match,
                                   bits_for_sequences,
                                   reinterpret_cast<uint8_t *>
                                     (encoded_match_list.data()),
                                   encoded_match_list.size(),
                                   reversed_byte_order);
    remove_duplicates<BytesUnitMatch>(&encoded_match_list);
  }
  size_t number_of_all_matches_get(void) const noexcept
  {
    return number_of_all_matches;
  }
  size_t number_of_seeds_get(void) const noexcept
  {
    return number_of_seeds;
  }
  int bits_for_sequences_get(void) const noexcept
  {
    return bits_for_sequences;
  }
  void statistics(void) const noexcept
  {
    std::cout << "# bits for sequences\t" << bits_for_sequences_get()
              << std::endl;
    std::cout << "# number of seeds\t" << number_of_seeds_get() << std::endl;
    std::cout << "# number of matches\t" << number_of_all_matches_get()
              << std::endl;
    std::cout << "# number of unique matches\t" << encoded_match_list.size()
              << std::endl;
  }
  /* This is currently not used */
  std::pair<size_t,size_t> gaps_of_adjacent(const BytesUnitMatch
                                               *segment_match_list,
                                            size_t i,size_t j) const noexcept
  {
    /* The following definition must be consistent with the definition in
       output_matches */
    constexpr const int ref_pos_idx = ref_idx == 0 ? 2 : 3;
    constexpr const int query_pos_idx = ref_idx == 0 ? 3 : 2;
    auto p = segment_match_list[i];
    auto m = segment_match_list[j];
    const uint64_t p_endpos0 = p.template decode_at<ref_pos_idx>(match_packer);
    const uint64_t p_endpos1
      = p.template decode_at<query_pos_idx>(match_packer);
    const uint64_t m_endpos0 = m.template decode_at<ref_pos_idx>(match_packer);
    const uint64_t m_endpos1
      = m.template decode_at<query_pos_idx>(match_packer);
    const uint64_t m_match_length
      = minimum_mem_length + m.template decode_at<4>(match_packer);
    const uint64_t m_startpos0 = m_endpos0 - m_match_length + 1,
                   m_startpos1 = m_endpos1 - m_match_length + 1;
    assert(p_endpos0 < m_startpos0 && p_endpos1 < m_startpos1);
    return std::make_pair(static_cast<size_t>(m_startpos0 - p_endpos0 - 1),
                          static_cast<size_t>(m_startpos1 - p_endpos1 - 1));
  }
  bool all_same_segment(void) const noexcept
  {
    return ref_multiseq->sequences_number_get() == 1 &&
           query_multiseq->sequences_number_get() == 1;
  }
  size_t size(void) const noexcept
  {
    return encoded_match_list.size();
  }
  bool same_segment(size_t i, size_t j) const noexcept
  {
    assert(i < size() && j < size());
    const BytesUnitMatch &a = encoded_match_list[i];
    const BytesUnitMatch &b = encoded_match_list[j];
    return (a.template decode_at<ref_idx>(match_packer) ==
            b.template decode_at<ref_idx>(match_packer)) &&
           (a.template decode_at<query_idx>(match_packer) ==
            b.template decode_at<query_idx>(match_packer));
  }
  size_t ref_seqnum_get(size_t this_idx) const noexcept
  {
    assert(this_idx < size());
    return static_cast<size_t>(encoded_match_list[this_idx]
                               .template decode_at<ref_idx>(match_packer));
  }
  size_t query_seqnum_get(size_t this_idx) const noexcept
  {
    assert(this_idx < size());
    return static_cast<size_t>(encoded_match_list[this_idx]
                               .template decode_at<query_idx>(match_packer));
  }
  size_t ref_endpos_get(size_t this_idx) const noexcept
  {
    assert(this_idx < size());
    return static_cast<size_t>(encoded_match_list[this_idx]
                               .template decode_at<ref_pos_idx>(match_packer));
  }
  size_t query_endpos_get(size_t this_idx) const noexcept
  {
    assert(this_idx < size());
    return static_cast<size_t>
                      (encoded_match_list[this_idx]
                       .template decode_at<query_pos_idx>(match_packer));
  }
  size_t length_get(size_t this_idx) const noexcept
  {
    assert(this_idx < size());
    return minimum_mem_length +
           static_cast<size_t>(encoded_match_list[this_idx]
                               .template decode_at<4>(match_packer));
  }
  size_t order_endpos_get(size_t this_idx) const noexcept
  {
    assert(this_idx < size());
    return static_cast<size_t>(encoded_match_list[this_idx]
                               .template decode_at<2>(match_packer));
  }
  std::pair<uint64_t,uint64_t> colinear_match_pair(size_t i,size_t j)
    const noexcept
  {
    auto p = encoded_match_list[i];
    auto m = encoded_match_list[j];
    const uint64_t j_endpos0
      = m.template decode_at<ref_pos_idx>(match_packer);
    const uint64_t j_endpos1
      = m.template decode_at<query_pos_idx>(match_packer);
    const size_t j_match_length = minimum_mem_length +
                                  m.template decode_at<4>(match_packer);
    assert(j_endpos0 + 1 >= j_match_length &&
           j_endpos1 + 1 >= j_match_length);
    const uint64_t j_startpos0 = j_endpos0 + 1 - j_match_length,
                   j_startpos1 = j_endpos1 + 1 - j_match_length;

    const uint64_t p_endpos0 = p.template decode_at<ref_pos_idx>(match_packer);

    if (p_endpos0 < j_startpos0)
    {
      const uint64_t p_endpos1
        = p.template decode_at<query_pos_idx>(match_packer);
      if (p_endpos1 < j_startpos1)
      {
        const uint64_t gap0 = j_startpos0 - p_endpos0 - 1,
                       gap1 = j_startpos1 - p_endpos1 - 1;
        assert(gap0 > 0 || gap1 > 0);
        return std::pair<uint64_t,uint64_t>(gap0,gap1);
      }
    }
    return std::pair<uint64_t,uint64_t>(0,0);
  }
};
#endif
