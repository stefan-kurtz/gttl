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
#include "utilities/ska_lsb_radix_sort.hpp"
#include "utilities/remove_duplicates.hpp"
#include "utilities/bytes_unit.hpp"
#include "sequences/gttl_multiseq.hpp"
#include "sequences/matching_characters.hpp"

template<bool check_bounds,bool (*matching_method)(char,char)>
static inline std::pair<size_t,size_t> maximize_on_both_ends(
                                               const char *ref_seq,
                                               size_t ref_seq_len,
                                               size_t ref_seq_abs_start_pos,
                                               const char *query_seq,
                                               size_t query_seq_len,
                                               size_t query_seq_abs_start_pos,
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

  assert(ref_seq_abs_start_pos + seed_size - 1 < ref_seq_len &&
         query_seq_abs_start_pos + seed_size - 1 < query_seq_len);
  const char *ref_seq_bwd, *query_seq_bwd;
  for (ref_seq_bwd = ref_seq + ref_seq_abs_start_pos - 1,
       query_seq_bwd = query_seq + query_seq_abs_start_pos - 1;
       pointers_leq(ref_seq,ref_seq_bwd) &&
       pointers_leq(query_seq,query_seq_bwd) &&
       matching_method(*ref_seq_bwd,*query_seq_bwd);
       ref_seq_bwd--, query_seq_bwd--)
       /* Nothing */ ;
  assert(ref_seq + ref_seq_abs_start_pos >= (ref_seq_bwd+1));
  const size_t left_extend
    = static_cast<size_t>(ref_seq + ref_seq_abs_start_pos - 1
                                  - (ref_seq_bwd+1) + 1);

  const char *ref_seq_fwd, *query_seq_fwd;
  for (ref_seq_fwd = ref_seq + ref_seq_abs_start_pos + seed_size,
       query_seq_fwd = query_seq + query_seq_abs_start_pos + seed_size;
       pointers_leq(ref_seq_fwd,ref_seq + ref_seq_len - 1) &&
       pointers_leq(query_seq_fwd,query_seq + query_seq_len - 1) &&
       matching_method(*ref_seq_fwd,*query_seq_fwd);
       ref_seq_fwd++, query_seq_fwd++)
         /* Nothing */;
  assert(ref_seq_fwd >= ref_seq + ref_seq_abs_start_pos + seed_size);
  const size_t right_extend
    = static_cast<size_t>(ref_seq_fwd - (ref_seq + ref_seq_abs_start_pos
                                                 + seed_size));
  return std::make_pair(left_extend,right_extend);
}

#define MAXIMIZE_ON_BOTH_ENDS(CHECK_BOUNDS,MATCHING_CHARACTERS)\
        std::tie(left_extend,right_extend) \
          = maximize_on_both_ends<CHECK_BOUNDS,MATCHING_CHARACTERS>\
                                 (ref_seq,\
                                  ref_seq_len,\
                                  pp.startpos0,\
                                  query_seq,\
                                  query_seq_len,\
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
    const uint64_t maximum_storable_match_length
      = gttl_bits2maxvalue<uint64_t>(remaining_bits_for_length);

    /* For the following loop to work the
       seed_enumerator must provide iterators begin() and end() which
       when dereferenced deliver an instance of
       struct SortedMatchListPositionPair, or some struct which is compatible.
    */
    for (auto const &pp : seed_enumerator)
    {
      if constexpr (seed_output)
      {
        std::cout << "# seed\t"
                  << pp.seqnum0 << "\t"
                  << pp.startpos0 << "\t"
                  << pp.seqnum1 << "\t"
                  << pp.startpos1 << std::endl;
      }
      const char *ref_seq = ref_multiseq->sequence_ptr_get(pp.seqnum0);
      const char *query_seq = query_multiseq->sequence_ptr_get(pp.seqnum1);
      const size_t ref_seq_len = ref_multiseq->sequence_length_get(pp.seqnum0);
      const size_t query_seq_len
        = query_multiseq->sequence_length_get(pp.seqnum1);
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
                                       (ref_seq + pp.startpos0,
                                        query_seq + pp.startpos1,
                                        qgram_length) == 0);
        } else
        {
          seed_okay = (count_mismatches<matching_characters>
                                       (ref_seq + pp.startpos0,
                                        query_seq + pp.startpos1,
                                        qgram_length) == 0);
        }
      }
      if (seed_okay)
      {
        assert(pp.startpos0 >= left_extend && pp.startpos1 >= left_extend);
        const size_t this_match_length = left_extend + qgram_length +
                                         right_extend;
        uint64_t length_stored;
        if (this_match_length >
            minimum_mem_length + maximum_storable_match_length)
        {
          length_stored = maximum_storable_match_length;
          if (sizeof_unit_match == 8)
          {
            StrFormat msg("cannot store match of length %lu in 8 bytes, resort "
                          "to using 9 bytes of space for each MEM",
                          this_match_length);
            throw msg.str();
          } else
          {
            StrFormat msg("cannot store match of length %lu in 9 bytes, please "
                          "inform the developer",
                          this_match_length);
            throw msg.str();
          }
        } else
        {
          assert(this_match_length >= minimum_mem_length);
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
  size_t size(void) const noexcept
  {
    return encoded_match_list.size();
  }
  size_t ref_seqnum_get(size_t idx) const noexcept
  {
    assert(idx < size());
    return static_cast<size_t>(encoded_match_list[idx]
                               .template decode_at<ref_idx>(match_packer));
  }
  size_t query_seqnum_get(size_t idx) const noexcept
  {
    assert(idx < size());
    return static_cast<size_t>(encoded_match_list[idx]
                               .template decode_at<query_idx>(match_packer));
  }
  size_t ref_endpos_get(size_t idx) const noexcept
  {
    assert(idx < size());
    return static_cast<size_t>(encoded_match_list[idx]
                               .template decode_at<ref_pos_idx>(match_packer));
  }
  size_t query_endpos_get(size_t idx) const noexcept
  {
    assert(idx < size());
    return static_cast<size_t>
                      (encoded_match_list[idx]
                       .template decode_at<query_pos_idx>(match_packer));
  }
  size_t length_get(size_t idx) const noexcept
  {
    assert(idx < size());
    return minimum_mem_length +
           static_cast<size_t>(encoded_match_list[idx]
                               .template decode_at<4>(match_packer));
  }
  size_t order_endpos_get(size_t idx) const noexcept
  {
    assert(idx < size());
    return static_cast<size_t>(encoded_match_list[idx]
                               .template decode_at<2>(match_packer));
  }
  bool same_segment(size_t i, size_t j) const noexcept
  {
    return (this->ref_seqnum_get(i) == this->ref_seqnum_get(j) &&
            this->query_seqnum_get(i) == this->query_seqnum_get(j));
  }
  std::pair<uint64_t,uint64_t> colinear_match_pair(size_t i,size_t j)
    const noexcept
  {
    const uint64_t ref_endpos_j = this->ref_endpos_get(j);
    const uint64_t query_endpos_j = this->query_endpos_get(j);
    const size_t j_match_length = this->length_get(j);
    assert(ref_endpos_j + 1 >= j_match_length &&
           query_endpos_j + 1 >= j_match_length);
    const uint64_t ref_startpos_j = ref_endpos_j + 1 - j_match_length,
                   query_startpos_j = query_endpos_j + 1 - j_match_length;

    const uint64_t ref_endpos_i = this->ref_endpos_get(i);
    if (ref_endpos_i < ref_startpos_j)
    {
      const uint64_t query_endpos_i = this->query_endpos_get(i);
      if (query_endpos_i < query_startpos_j)
      {
        const uint64_t ref_gap = ref_startpos_j - ref_endpos_i - 1,
                       query_gap = query_startpos_j - query_endpos_i - 1;
        assert(ref_gap > 0 || query_gap > 0);
        return std::make_pair(ref_gap,query_gap);
      }
    }
    return std::make_pair(uint64_t(0),uint64_t(0));
  }
  void statistics(FILE *out_fp) const noexcept
  {
    fprintf(out_fp,"# bits for sequences\t%d\n",bits_for_sequences);
    fprintf(out_fp,"# number of seeds\t%lu\n",number_of_seeds_get());
    fprintf(out_fp,"# number of all matches\t%lu\n",
            number_of_all_matches_get());
    fprintf(out_fp,"# number of unique matches\t%lu\n",size());
  }
  bool all_same_segment(void) const noexcept
  {
    return ref_multiseq->sequences_number_get() == 1 &&
           query_multiseq->sequences_number_get() == 1;
  }
  /* This is currently not used */
  std::pair<size_t,size_t> gaps_of_adjacent(const BytesUnitMatch
                                               *segment_match_list,
                                            size_t i,size_t j) const noexcept
  {
    /* The following definition must be consistent with the definition in
       output_matches */
    const BytesUnitMatch &p = segment_match_list[i];
    const BytesUnitMatch &m = segment_match_list[j];
    const uint64_t ref_endpos_p
      = p.template decode_at<ref_pos_idx>(match_packer);
    const uint64_t query_endpos_p
      = p.template decode_at<query_pos_idx>(match_packer);
    const uint64_t ref_endpos_m
      = m.template decode_at<ref_pos_idx>(match_packer);
    const uint64_t query_endpos_m
      = m.template decode_at<query_pos_idx>(match_packer);
    const uint64_t m_match_length
      = minimum_mem_length + m.template decode_at<4>(match_packer);
    const uint64_t ref_startpos_m = ref_endpos_m - m_match_length + 1,
                   query_startpos_m = query_endpos_m - m_match_length + 1;
    assert(ref_endpos_p < ref_startpos_m && query_endpos_p < query_startpos_m);
    return std::make_pair(static_cast<size_t>(ref_startpos_m - ref_endpos_p
                                                             - 1),
                          static_cast<size_t>(query_startpos_m - query_endpos_p
                                                               - 1));
  }
};
#endif
