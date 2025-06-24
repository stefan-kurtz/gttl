#ifndef CLOSE_GAPS_HPP
#define CLOSE_GAPS_HPP
#include <cstddef>
#include <cstdint>
#include <cassert>
#include "sequences/matching_characters.hpp"
#include "sequences/outsenseedist_unr.hpp"

template<typename CharType,bool (*match_method)(CharType,CharType)>
static inline size_t single_symbol_edist(CharType single_cc,
                                         const CharType *other,
                                         size_t gap_len)
{
  size_t matches = 0;
  for (size_t idx = 0; idx < gap_len; idx++)
  {
    matches += match_method(single_cc,other[idx]);
  }
  return (matches > 0 ? gap_len - 1 : gap_len);
}

template<typename CharType,bool (*match_method)(CharType,CharType)>
static inline size_t equal_length_sequences_edist(const CharType *useq,
                                                  const CharType *vseq,
                                                  size_t gap_len,
                                                  size_t useqnum,
                                                  size_t vseqnum)
{
  assert(gap_len > 2);
  /* As we know tha the gap is of length > 2 and the boundaries
     of the gap are distinct, we ony determine the number of mismatches
     of the inner part of the gap. The number of mismatches gives an upper
     bound on the edit distance, which we use to determin the d_max parameter.
  */
  const size_t mismatches = count_mismatches<CharType,match_method>
                                            (useq + 1,
                                             vseq + 1,
                                             gap_len - 2);
  const size_t d_max = mismatches + 2;
  const size_t edist
    = fastedist_unrolled_same_seq_length<CharType,size_t,
                                         lcplen_fwd<match_method,false>>
                                        (d_max,
                                         useq,
                                         gap_len,
                                         vseq,
                                         gap_len,
                                         useqnum,
                                         vseqnum);
  assert(edist <= d_max);
  return edist;
}

template<typename CharType,bool (*match_method)(CharType,CharType)>
static inline size_t different_length_sequences_edist(
                                               const CharType *useq,
                                               const CharType *vseq,
                                               size_t ulen,
                                               size_t vlen,
                                               size_t useqnum,
                                               size_t vseqnum)
{
  if (ulen < vlen)
  {
    return fastedist_unrolled<CharType,size_t,lcplen_fwd<match_method,false>>
                             (vlen,
                              useq,
                              ulen,
                              vseq,
                              vlen,
                              useqnum,
                              vseqnum);
  }
  return fastedist_unrolled<CharType,size_t,lcplen_fwd<match_method,false>>
                           (ulen,
                            vseq,
                            vlen,
                            useq,
                            ulen,
                            useqnum,
                            vseqnum);
}

#ifndef NDEBUG
template<typename CharType,bool (*match_method)(CharType,CharType)>
static void check_gap_boundaries(const CharType *useq,
                                 const CharType *vseq,
                                 size_t p_uendpos,
                                 size_t ulen,
                                 size_t p_vendpos,
                                 size_t vlen)
{
  assert(!match_method(useq[p_uendpos + 1],
                       vseq[p_vendpos + 1]) &&
         !match_method(useq[p_uendpos + ulen],
                       vseq[p_vendpos + vlen]));
}
#endif

template<class SeedTable,
         typename CharType,
         bool (*match_method)(CharType,CharType)>
static size_t all_sequence_types_edist(const SeedTable &seed_table,
                                       size_t previous_idx,
                                       size_t ulen,
                                       size_t vlen,
                                       const CharType *useq,
                                       const CharType *vseq,
                                       size_t useqnum,
                                       size_t vseqnum)
{
  if (ulen == 0)
  {
    return vlen;
  }
  if (vlen == 0)
  {
    return ulen;
  }
  if (ulen == vlen)
  {
    if (ulen <= 2)
    {
      return ulen;
    }
    const uint64_t p_uendpos = seed_table.ref_endpos_get(previous_idx);
    const uint64_t p_vendpos = seed_table.query_endpos_get(previous_idx);
#ifndef NDEBUG
    check_gap_boundaries<CharType,match_method>
                        (useq,
                         vseq,
                         p_uendpos,
                         ulen,
                         p_vendpos,
                         vlen);
#endif
    const size_t edist
      = equal_length_sequences_edist<CharType,match_method>
                                    (useq + p_uendpos + 1,
                                     vseq + p_vendpos + 1,
                                     ulen,
                                     useqnum,
                                     vseqnum);
    assert(edist >= 2);
    return edist;
  }
  const uint64_t p_uendpos = seed_table.ref_endpos_get(previous_idx);
  const uint64_t p_vendpos = seed_table.query_endpos_get(previous_idx);
#ifndef NDEBUG
  check_gap_boundaries<CharType,match_method>
                      (useq,
                       vseq,
                       p_uendpos,
                       ulen,
                       p_vendpos,
                       vlen);
#endif
  if (ulen == 1)
  {
    return single_symbol_edist<CharType,match_method>
                              (useq[p_uendpos + 1],
                               vseq + p_vendpos + 1,
                               vlen);
  }
  if (vlen == 1)
  {
    return single_symbol_edist<CharType,match_method>
                              (vseq[p_vendpos + 1],
                               useq + p_uendpos + 1,
                               ulen);
  }
  const size_t edist
    = different_length_sequences_edist<CharType,match_method>
                                      (useq + p_uendpos + 1,
                                       vseq + p_vendpos + 1,
                                       ulen,
                                       vlen,
                                       useqnum,
                                       vseqnum);
  assert(edist >= 2);
  return edist;
}

template<class SeedTable,
         class ChainerClass,
         typename CharType,
         bool (*match_method)(CharType,CharType)>
static size_t close_gaps_between_chain_elements(const SeedTable &seed_table,
                                                size_t segment_start,
                                                const ChainerClass
                                                  &local_chainer,
                                                const CharType *useq,
                                                const CharType *vseq,
                                                size_t from_element,
                                                size_t chain_len,
                                                size_t useqnum,
                                                size_t vseqnum)
{
  if (chain_len == 1)
  {
    return 0;
  }
  size_t j = from_element,
         distance = 0;
  for (size_t current_length = 1; current_length < chain_len; current_length++)
  {
    const size_t ulen = local_chainer.ref_gap_length_get(j),
                 vlen = local_chainer.query_gap_length_get(j);
    j = static_cast<size_t>(local_chainer.predecessor_get(j));
    distance += all_sequence_types_edist<SeedTable,CharType,match_method>
                                        (seed_table,
                                         segment_start + j,
                                         ulen,
                                         vlen,
                                         useq,
                                         vseq,
                                         useqnum,
                                         vseqnum);
  }
  return distance;
}
#endif
