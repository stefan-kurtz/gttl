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
#ifndef FILL_BU_SUFTAB_HPP
#define FILL_BU_SUFTAB_HPP
#include <cstddef>
#include <cassert>
#include <string>
#include <cstdint>
#include <vector>
#include <cstdlib>
#include "utilities/read_vector.hpp"
#include "utilities/ordered_integer_sequence.hpp"
#include "utilities/gttl_binary_write.hpp"
#include "utilities/memory_tracker.hpp"
#include "sequences/gttl_multiseq.hpp"
#include "inverse_suftab_iter.hpp"

template <typename intset_basetype>
static OrderedIntegerSequence<intset_basetype> separator_positions_get(
                                                   const GttlMultiseq *multiseq)
{
  const size_t sequences_number = multiseq->sequences_number_get();
  const size_t totallength
    = multiseq->sequences_total_length_get() + sequences_number - 1;

  OrderedIntegerSequence<intset_basetype>
    separator_positions(totallength, sequences_number);
  size_t sum_length = multiseq->sequence_length_get(0);
  separator_positions.append(sum_length);
  for (size_t idx = 1; idx < sequences_number; idx++)
  {
    const size_t length = multiseq->sequence_length_get(idx);
    sum_length += (length + 1);
    separator_positions.append(sum_length);
  }
  assert(sum_length == totallength);
  return separator_positions;
}

template <typename SuftabBaseType,typename T_bp, typename T_bu,
          typename intset_basetype>
static void fill_bu_suftab_MULTISEQ(const SuftabBaseType *suftab,
                                    const T_bp &bp,
                                    const GttlMultiseq *multiseq,
                                    size_t totallength,
                                    bool verbose,
                                    const std::string &indexname)
{
  const OrderedIntegerSequence<intset_basetype> separator_positions =
                               separator_positions_get<intset_basetype>(
                                                            multiseq);

  const std::string bsf_filename(indexname + ".bsf");
  BinaryFileWriter<T_bu> bw(bsf_filename);
  if (verbose)
  {
    separator_positions.show_info();
  }
  for (size_t idx = 0; idx <= totallength; idx++)
  {
    const size_t seqnum = separator_positions.pos2seqnum(suftab[idx]);
    size_t pos;
    if (seqnum > 0)
    {
      const size_t sectionnum_upperbound
        = separator_positions.section_number_get(suftab[idx]);
      const size_t sectionnum
        = separator_positions
            .backward_sec_idx_largest_leq(seqnum - 1,
                                          sectionnum_upperbound);
      const size_t sequence_start
        = separator_positions.get_element_at(sectionnum, seqnum - 1) + 1;
      assert(suftab[idx] >= sequence_start);
      pos = suftab[idx] - sequence_start;
    } else
    {
      pos = suftab[idx];
    }
    bw.append(T_bu(bp, {seqnum, pos}));
  }
}

template <typename SuftabBaseType,typename T_bp, typename T_bu>
static void fill_bu_suftab_MULTISEQ_linear(GttlMemoryTracker *memory_tracker,
                                           const SuftabBaseType *suftab,
                                           T_bu *bu_suftab,
                                           const T_bp &bp,
                                           const GttlMultiseq *multiseq,
                                           size_t totallength)
{
  SuftabBaseType *const inverse_suftab = new SuftabBaseType[totallength + 1];
  memory_tracker->track(inverse_suftab,__FILE__,__LINE__,
                        (totallength + 1) * sizeof(SuftabBaseType));
  for (size_t idx = 0; idx <= totallength; idx++)
  {
    inverse_suftab[suftab[idx]] = idx;
  }
  const size_t sequences_number = multiseq->sequences_number_get();
  size_t idx = 0;
  for (size_t seqnum = 0; seqnum < sequences_number; seqnum++)
  {
    for (size_t relpos = 0; relpos < multiseq->sequence_length_get(seqnum) + 1;
         relpos++)
    {
      bu_suftab[inverse_suftab[idx++]] = T_bu(bp, {seqnum, relpos});
    }
  }
  memory_tracker->untrack(inverse_suftab,__FILE__,__LINE__);
  delete[] inverse_suftab;
}

template <typename SuftabBaseType,typename T_bp, typename T_bu>
static T_bu *fill_bu_suftab_MULTISEQ_linear(GttlMemoryTracker *memory_tracker,
                                            const std::string &indexname,
                                            const T_bp &bp,
                                            const GttlMultiseq *multiseq,
                                            size_t totallength)
{
  const InverseSuftabReader<SuftabBaseType> inverse_suftab_reader(
                               memory_tracker, indexname, totallength);
  const size_t sequences_number = multiseq->sequences_number_get();
  auto inverse_suftab_iter = inverse_suftab_reader.begin();
  T_bu *const bu_suftab         = new T_bu[totallength + 1];
  memory_tracker->track(bu_suftab,__FILE__,__LINE__,(totallength + 1)
                                                    * sizeof (T_bu));
  for (size_t seqnum = 0; seqnum < sequences_number; seqnum++)
  {
    for (size_t relpos = 0; relpos < multiseq->sequence_length_get(seqnum) + 1;
         relpos++)
    {
      bu_suftab[*inverse_suftab_iter] = T_bu(bp, {seqnum, relpos});
      ++inverse_suftab_iter;
    }
  }
  return bu_suftab;
}

template <typename SuftabBaseType,typename T_bu, typename T_bp>
static void check_bu_suftab_MULTISEQ(const SuftabBaseType *suftab,
                                     const GttlMultiseq *multiseq,
                                     const T_bp &bp,
                                     size_t totallength,
                                     const std::string &filename)
{
  const uint8_t padding_char = multiseq->padding_char_get();
  const OrderedIntegerSequence<uint32_t> separator_positions_uint =
                               separator_positions_get<uint32_t>(multiseq);
  std::vector<T_bu> bu_suftab = gttl_read_vector<T_bu>(filename);

  for (size_t idx = 0; idx <= totallength; idx++)
  {
    if (multiseq->sequence_char_get(suftab[idx]) != padding_char)
    {
      const size_t seqnum = bu_suftab[idx].template decode_at<0>(bp);
      const size_t rel_pos = bu_suftab[idx].template decode_at<1>(bp);
      size_t abs_pos;
      if (seqnum > 0)
      {
        const size_t sequence_start
          = separator_positions_uint.get_element_at(seqnum - 1) + 1;
        abs_pos = sequence_start + rel_pos;
      } else
      {
        abs_pos = rel_pos;
      }
      if (suftab[idx] != abs_pos)
      {
        fprintf(stderr,"%s: suftab[%zu] = %zu != %zu derived from "
                       "bu_suftab[%zu]",
                     __func__,idx,static_cast<size_t>(suftab[idx]),abs_pos,idx);
        exit(EXIT_FAILURE);
      }
    }
  }
}

template <typename SuftabBaseType,typename T_bp, typename T_bu>
static void fill_bu_suftab_PLAINSEQ(const SuftabBaseType *suftab,
                                    const T_bp &bp,
                                    size_t totallength,
                                    const std::string &indexname)
{
  const std::string bsf_filename(indexname + ".bsf");
  BinaryFileWriter<T_bu> bw(bsf_filename);
  for (size_t idx = 0; idx <= totallength; idx++)
  {
    bw.append(T_bu(bp, {0, suftab[idx]}));
  }
}
#endif
