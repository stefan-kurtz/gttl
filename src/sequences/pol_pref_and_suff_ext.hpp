#ifndef POL_PREF_AND_SUFF_EXT_HPP
#define POL_PREF_AND_SUFF_EXT_HPP
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include "sequences/polished_points.hpp"
#include "sequences/align_polish.hpp"
#include "sequences/pol_pref_or_suff_ext.hpp"

template<class MatchClass>
static bool combine_suffix_prefix_extensions(
                               MatchClass *extended_match,
                               const MatchClass &inner_match,
                               size_t aligned_len_threshold,
                               double maximum_error_percentage,
                               const PolishedPoints &pps_suffix,
                               const PolishedPoints &pps_prefix)
{
  const size_t upper_bound_aligned_len
    = inner_match.aligned_len_get() + pps_suffix.longest_aligned_len_get()
                                    + pps_prefix.longest_aligned_len_get();
#undef SKDEBUG
#ifdef SKDEBUG
  printf("upper_bound_aligned_len=%zu=%zu+%zu+%zu\n",
         upper_bound_aligned_len,
         inner_match.aligned_len_get(),
         pps_suffix.longest_aligned_len_get(),
         pps_prefix.longest_aligned_len_get());
#endif
  if (upper_bound_aligned_len >= aligned_len_threshold)
  {
    const size_t lower_bound_distance
      = inner_match.distance_get() +
        pps_suffix.smallest_unit_cost_get() +
        pps_prefix.smallest_unit_cost_get();
#ifdef SKDEBUG
    printf("lower_bound_distance=%zu=%zu+%zu+%zu\n",
            lower_bound_distance,
            inner_match.distance_get(),
            pps_suffix.smallest_unit_cost_get(),
            pps_prefix.smallest_unit_cost_get());
    printf("error_percentage=%.2f\n",
           error_percentage_get(lower_bound_distance,upper_bound_aligned_len));
#endif
    if (error_percentage_get(lower_bound_distance,upper_bound_aligned_len) <=
        maximum_error_percentage)
    {
      for (auto &&pp_s : pps_suffix)
      {
        for (auto &&pp_p : pps_prefix)
        {
          const size_t this_aligned_len = inner_match.aligned_len_get() +
                                          pp_s.aligned_len_get() +
                                          pp_p.aligned_len_get();
          if (this_aligned_len >= aligned_len_threshold)
          {
            const size_t this_distance
              = inner_match.distance_get() + pp_s.distance_get()
                                           + pp_p.distance_get();
            if (error_percentage_get(this_distance,this_aligned_len)
                <= maximum_error_percentage)
            {
#ifdef SKDEBUG
              printf("s_endpos=%zu+%zu\n",inner_match.s_endpos_get(),
                                          pp_p.row_get());
              printf("q_endpos=%zu+%zu\n",inner_match.q_endpos_get(),
                                          pp_p.column_get());
              printf("s_len=%zu+%zu+%zu\n",pp_s.row_get(),
                                           inner_match.s_len_get(),
                                           pp_p.row_get());
              printf("q_len=%zu+%zu+%zu\n",pp_s.column_get(),
                                           inner_match.q_len_get(),
                                           pp_p.column_get());
#endif
              extended_match->s_endpos_set(inner_match.s_endpos_get() +
                                           pp_p.row_get());
              extended_match->q_endpos_set(inner_match.q_endpos_get() +
                                           pp_p.column_get());
              extended_match->s_len_set(inner_match.s_len_get() +
                                        pp_s.row_get() + pp_p.row_get());
              extended_match->q_len_set(inner_match.q_len_get() +
                                        pp_s.column_get() + pp_p.column_get());
              extended_match->distance_set(this_distance);
              return true;
            }
          }
        }
      }
    }
  }
  if (inner_match.aligned_len_get() >= aligned_len_threshold &&
      inner_match.distance_get() <= maximum_error_percentage)
  {
    extended_match->s_endpos_set(inner_match.s_endpos_get());
    extended_match->q_endpos_set(inner_match.q_endpos_get());
    extended_match->s_len_set(inner_match.s_len_get());
    extended_match->q_len_set(inner_match.q_len_get());
    extended_match->distance_set(inner_match.distance_get());
    return true;
  }
  return false;
}

template<class MatchClass>
static bool polished_prefix_and_suffix_extension(
                                  MatchClass *extended_match,
                                  const AlignmentPolishing &alignment_polishing,
                                  size_t minimum_als_length,
                                  double maximum_error_percentage,
                                  const char *ref_seq,
                                  size_t ref_len,
                                  const char *query_seq,
                                  size_t query_len,
                                  size_t ref_seqnum,
                                  size_t query_seqnum,
                                  const MatchClass &inner_match)
{
  assert(ref_len > inner_match.s_endpos_get() &&
         query_len > inner_match.q_endpos_get());

  PolishedPoints best_polished_points_suffix(inner_match.distance_get(),
                                             inner_match.aligned_len_get());
  prefix_or_suffix_extension<lcslen_bwd<matching_characters_wc>>
                            (&best_polished_points_suffix,
                             alignment_polishing,
                             ref_seq,
                             inner_match.s_startpos_get(), /* length */
                             query_seq,
                             inner_match.q_startpos_get(), /* length */
                             ref_seqnum,
                             query_seqnum);
  PolishedPoints best_polished_points_prefix(inner_match.distance_get(),
                                             inner_match.aligned_len_get());
  prefix_or_suffix_extension<lcplen_fwd<matching_characters_wc,false>>
                            (&best_polished_points_prefix,
                             alignment_polishing,
                             ref_seq + inner_match.s_endpos_get() + 1,
                             ref_len - 1 - inner_match.s_endpos_get(),
                             query_seq + inner_match.q_endpos_get() + 1,
                             query_len - 1 - inner_match.q_endpos_get(),
                             ref_seqnum,
                             query_seqnum);

  const size_t aligned_len_threshold = 2 * minimum_als_length;
#ifdef SKDEBUG
  printf("best_polished_points_suffix\t%s\n",
         best_polished_points_suffix.to_string().c_str());
  printf("best_polished_points_prefix\t%s\n",
         best_polished_points_prefix.to_string().c_str());
  printf("inner_match.distance\t%zu\n",inner_match.distance_get());
#endif
  return combine_suffix_prefix_extensions<MatchClass>
                                         (extended_match,
                                          inner_match,
                                          aligned_len_threshold,
                                          maximum_error_percentage,
                                          best_polished_points_suffix,
                                          best_polished_points_prefix);
}
#endif
