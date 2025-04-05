#ifndef SW_COMPARATOR_HPP
#define SW_COMPARATOR_HPP

#include <cstdint>
#include <cassert>
#include "alignment/loc_align_coords.hpp"
#include "alignment/ssw.hpp"

template<class SWProcessPair,class SWProcessResultShared,
         class SWProcessResultThreadRelated>
class SWcomparator
{
  SWProcessResultThreadRelated *sw_process_result_thread_related;
  size_t alphasize;
  bool try_reverse_strand;
  const int8_t *const *scorematrix2D;
  int8_t smallest_score;
  int8_t gap_open_penalty;
  int8_t gap_extension_penalty;
  SSWprofile *ssw_profile;
  SSWresources ssw_resources;
  size_t current_i;
  bool compute_only_end;
  const SWProcessPair &process_this_pair;
  const SWProcessResultShared &sw_process_result_shared;

  public:
  SWcomparator(SWProcessResultThreadRelated *_sw_process_result_thread_related,
               size_t _alphasize,
               bool _try_reverse_strand,
               const int8_t *const *_scorematrix2D,
               int8_t _smallest_score,
               int8_t _gap_open_penalty,
               int8_t _gap_extension_penalty,
               size_t max_seq_len,
               bool _compute_only_end,
               const SWProcessPair &_process_this_pair,
               const SWProcessResultShared &_sw_process_result_shared)
    : sw_process_result_thread_related(_sw_process_result_thread_related)
    , alphasize(_alphasize)
    , try_reverse_strand(_try_reverse_strand)
    , scorematrix2D(_scorematrix2D)
    , smallest_score(_smallest_score)
    , gap_open_penalty(_gap_open_penalty)
    , gap_extension_penalty(_gap_extension_penalty)
    , ssw_profile(nullptr)
    , ssw_resources(8 + 16,max_seq_len)
    , current_i(~size_t(0))
    , compute_only_end(_compute_only_end)
    , process_this_pair(_process_this_pair)
    , sw_process_result_shared(_sw_process_result_shared)
  { }
  ~SWcomparator(void)
  {
    delete ssw_profile;
  }
  void preprocess(size_t i,const std::string_view &db_seq)
  {
    if (ssw_profile != nullptr)
    {
      delete ssw_profile;
    }
    current_i = i;
    ssw_profile = new SSWprofile (alphasize,
                                  scorematrix2D[0],
                                  smallest_score,
                                  reinterpret_cast<const uint8_t *>
                                    (db_seq.data()),
                                  db_seq.size());
  }
  bool compare(size_t i,size_t j,const std::string_view &query_seq)
  {
    if (not process_this_pair.check(i,j))
    {
      return false;
    }
    assert(i == current_i);
    LocalAlignmentCoordinates la_coords;
    static constexpr const bool forward_strand = true;
    ssw_resources.reset8(ssw_len2segment_len<uint8_t>(query_seq.size()));
    ssw_resources.reset16(ssw_len2segment_len<uint16_t>(query_seq.size()));
    ssw_resources.reset32(ssw_len2segment_len<uint32_t>(query_seq.size()));
    std::tie(la_coords.raw_score,
             la_coords.ustart,
             la_coords.usubstringlength,
             la_coords.vstart,
             la_coords.vsubstringlength)
      = ssw_align<forward_strand>
                 (*ssw_profile,
                  &ssw_resources,
                  reinterpret_cast<const uint8_t *>(query_seq.data()),
                  query_seq.size(),
                  /* as ssw does not count gap_extension for
                  first symbol in gap, we add it to the
                  gap_open_penalty here */
                  static_cast<uint8_t>(gap_open_penalty +
                                       gap_extension_penalty),
                  static_cast<uint8_t>(gap_extension_penalty),
                  compute_only_end);
    la_coords.forward_strand = true;

    LocalAlignmentCoordinates la_coords_rc;
    if (try_reverse_strand)
    {
      std::tie(la_coords_rc.raw_score,
               la_coords_rc.ustart,
               la_coords_rc.usubstringlength,
               la_coords_rc.vstart,
               la_coords_rc.vsubstringlength)
        = ssw_align<not forward_strand>
                   (*ssw_profile,
                    &ssw_resources,
                    reinterpret_cast<const uint8_t *>(query_seq.data()),
                    query_seq.size(),
                    /* as ssw does not count gap_extension for
                    first symbol in gap, we add it to the
                    gap_open_penalty here */
                    static_cast<uint8_t>(gap_open_penalty +
                                         gap_extension_penalty),
                    static_cast<uint8_t>(gap_extension_penalty),
                    compute_only_end);
      la_coords_rc.forward_strand = false;
      if (la_coords_rc > la_coords)
      {
        return sw_process_result_shared.process(
                                              sw_process_result_thread_related,
                                              la_coords_rc,i,j);
      }
    }
    return sw_process_result_shared.process(sw_process_result_thread_related,
                                            la_coords,i,j);
  }
};
#endif
