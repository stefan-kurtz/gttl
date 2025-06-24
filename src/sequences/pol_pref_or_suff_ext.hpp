#ifndef POL_PREF_OR_SUFF_EXT_HPP
#define POL_PREF_OR_SUFF_EXT_HPP
#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <cassert>
#include "sequences/align_polish.hpp"
#include "sequences/lcs_lcp_len_type.hpp"
#include "sequences/prefix_suffix_extension.hpp"
#include "sequences/polished_points.hpp"
#include "sequences/front_value_match_history.hpp"

template<LcsLcpLenType suffix_or_prefix_match_len>
static void prefix_or_suffix_extension(PolishedPoints *best_polished_points,
                                       const AlignmentPolishing
                                         &alignment_polishing,
                                       const char *useq,
                                       size_t ulen,
                                       const char *vseq,
                                       size_t vlen,
                                       size_t useqnum,
                                       size_t vseqnum)
{
  static constexpr const bool track_eop = false;
  TrackPolishedPoints<FrontValueMatchHistory>
    track_polished_points(best_polished_points,alignment_polishing, ulen,vlen);
  prefix_or_suffix_extension_generic<track_eop,
                                     TrackPolishedPoints
                                       <FrontValueMatchHistory>,
                                     FrontValueMatchHistory,
                                     suffix_or_prefix_match_len>
                                     (&track_polished_points,
                                      useq,
                                      ulen,
                                      vseq,
                                      vlen,
                                      useqnum,
                                      vseqnum);
}
#endif
