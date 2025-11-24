#ifndef OUTSENSEEDIST_INPLACE_HPP
#define OUTSENSEEDIST_INPLACE_HPP
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <utility>
#include "sequences/lcs_lcp_len_type.hpp"

template<bool track_eop,class FrontValue,
         LcsLcpLenType suffix_or_prefix_match_len>
static inline void outsense_next_front_after_second_inplace(
                                                FrontValue *front,
                                                size_t d,
                                                const char *useq,
                                                size_t ulen,
                                                const char *vseq,
                                                size_t vlen,
                                                size_t seqnum0,
                                                size_t seqnum1)
{
  assert(d > 0);
  FrontValue insertion_value = front[0]; /* from previous diag -(d-1)
                                            => -d => DELETION */
  FrontValue bestfront = insertion_value;
  bestfront += 1; /* only increment row */
  front[0] = bestfront;
  if constexpr (track_eop)
  {
    front[0].deletion_set();
  }
  if (bestfront < ulen && bestfront < vlen + d)
  {
    assert(bestfront + 0 >= d);
    const size_t l = suffix_or_prefix_match_len(useq,bestfront + 0,
                                                vseq,(bestfront + 0) - d,
                                                ulen,vlen,
                                                seqnum0,seqnum1);
    front[0] += l; /* add match of length l */
  } else
  {
    front[0] += size_t(0);
  }
  FrontValue replacement_value = front[1];
  if (bestfront <= replacement_value)
  {
    bestfront = replacement_value;
    if constexpr (track_eop)
    {
      bestfront.deletion_set();
    }
    bestfront += 1; /* only increment row */
  } else
  {
    if constexpr (track_eop)
    {
      bestfront.mismatch_set();
      if (bestfront == replacement_value + 1)
      {
        bestfront.deletion_add();
      }
    }
  }
  front[1] = bestfront;
  if (bestfront < ulen && bestfront + 1 < vlen + d)
  {
    assert(bestfront + 1 >= d);
    const size_t l = suffix_or_prefix_match_len(useq,bestfront + 0,
                                                vseq,(bestfront + 1) - d,
                                                ulen,vlen,
                                                seqnum0,seqnum1);
    front[1] += l;
  } else
  {
    front[1] += size_t(0);
  }
  const size_t frontmaxidx = 2 * d;
  for (size_t idx=size_t(2); idx <= frontmaxidx; idx++)
  {
    bestfront = insertion_value;
    if constexpr (track_eop)
    {
      bestfront.insertion_set();
    }
    bestfront += 0; /* do not delete this line */
    if (idx <= frontmaxidx - 1)
    {
      if (bestfront <= replacement_value)
      {
        bestfront = replacement_value;
        if constexpr (track_eop)
        {
          bestfront.mismatch_set();
        }
        bestfront += 1; /* only increment row */
      } else
      {
        if constexpr (track_eop)
        {
          if (bestfront == replacement_value + 1)
          {
            bestfront.mismatch_add();
          }
        }
      }
    }
    if (idx <= frontmaxidx - 2)
    {
      if (bestfront <= front[idx])
      {
        bestfront = front[idx];
        if constexpr (track_eop)
        {
          bestfront.deletion_set();
        }
        bestfront += 1; /* only increment row */
      } else
      {
        if constexpr (track_eop)
        {
          if (bestfront == front[idx] + 1)
          {
            bestfront.deletion_add();
          }
        }
      }
    }
    if (idx < frontmaxidx)
    {
      insertion_value = replacement_value;
      replacement_value = front[idx];
    }
    front[idx] = bestfront;
    if (bestfront < ulen && bestfront + idx < vlen + d)
    {
      assert(bestfront + idx >= d);
      const size_t l = suffix_or_prefix_match_len(useq,bestfront + 0,
                                                  vseq,(bestfront + idx) - d,
                                                  ulen,vlen,
                                                  seqnum0,seqnum1);
      front[idx] += l;
    } else
    {
      front[idx] += size_t(0);
    }
  }
}

template<bool track_eop,class FrontValue,
         LcsLcpLenType suffix_or_prefix_match_len>
static inline void outsense_second_front_inplace(FrontValue *front,
                                                 const char *useq,
                                                 size_t ulen,
                                                 const char *vseq,
                                                 size_t vlen,
                                                 size_t seqnum0,
                                                 size_t seqnum1)
{
  front[1] = front[2] = front[0];
  if constexpr (track_eop)
  {
    front[0].deletion_set();
  }
  front[0] += 1; /* only increment row */
  if (front[0] < ulen && front[0] < vlen + 1)
  {
    assert((front[0] + 0) > 0);
    const size_t l = suffix_or_prefix_match_len(useq,front[0] + 0,
                                                vseq,(front[0] + 0) - 1,
                                                ulen,vlen,
                                                seqnum0,seqnum1);
    front[0] += l; /* increment row and set match_length */
  } else
  {
    front[0] += size_t(0);
  }
  if constexpr (track_eop)
  {
    front[1].mismatch_set();
  }
  front[1] += 1; /* only increment row */
  if (front[1] < ulen && front[1] < vlen)
  {
    const size_t l = suffix_or_prefix_match_len(useq,front[1] + 0,
                                                vseq,front[1] + 0,
                                                ulen,vlen,
                                                seqnum0,seqnum1);
    front[1] += l;
  } else
  {
    front[1] += size_t(0);
  }
  if constexpr (track_eop)
  {
    front[2].insertion_set();
  }
  front[2] += 0;
  if (front[2] < ulen && front[2] + 1 < vlen)
  {
    const size_t l = suffix_or_prefix_match_len(useq,front[2] + 0,
                                                vseq,front[2] + 1,
                                                ulen,vlen,
                                                seqnum0,seqnum1);
    front[2] += l;
  } else
  {
    front[2] += size_t(0);
  }
}

template<bool track_eop,class FrontValue,
         LcsLcpLenType suffix_or_prefix_match_len>
static inline void outsense_next_front_inplace(FrontValue *front,
                                               size_t d,
                                               const char *useq,
                                               size_t ulen,
                                               const char *vseq,
                                               size_t vlen,
                                               size_t seqnum0,
                                               size_t seqnum1)
{
  if (d == 1)
  {
    outsense_second_front_inplace<track_eop,FrontValue,
                                  suffix_or_prefix_match_len>
                                 (front,useq,ulen,vseq,vlen,seqnum0,seqnum1);
  } else
  {
    assert(d > 1);
    outsense_next_front_after_second_inplace<track_eop,FrontValue,
                                             suffix_or_prefix_match_len>
                                             (front,d,useq,ulen,vseq,vlen,
                                              seqnum0,seqnum1);
  }
}

template<typename char_type,class FrontValue,int previous_d,bool d_max_defined,
         LcsLcpLenType suffix_or_prefix_match_len>
static size_t fastedist_inplace_continue(const FrontValue *previousfront,
                                         size_t d_max,
                                         const char_type *useq,
                                         size_t ulen,
                                         const char_type *vseq,
                                         size_t vlen,
                                         size_t seqnum0,
                                         size_t seqnum1)
{
  assert(ulen > 0 && vlen > 0);
  size_t d;
  size_t allocated = 0;

  if constexpr (d_max_defined)
  {
    assert(d_max > 0 && std::cmp_less(previous_d, d_max));
    allocated = 2 * d_max + 1;
  } else
  {
    assert(d_max == 0);
    allocated = 2 * std::max(32,previous_d) + 1;
  }
  FrontValue *front
    = static_cast<FrontValue *>(malloc(allocated * sizeof *front));
  assert(front != nullptr);
  memcpy(front,previousfront,
         static_cast<size_t>(2 * previous_d + 1) * sizeof *front);
  for (d = static_cast<size_t>(previous_d)+1; /* Nothing */; d++)
  {
    if constexpr (!d_max_defined)
    {
      if (2 * d >= allocated)
      {
        allocated = std::max(2 * d + 1,(allocated * 12)/10 + size_t(128));
        front = static_cast<FrontValue *>
                           (realloc(front,allocated * sizeof *front));
        assert(front != nullptr);
      }
    }
    outsense_next_front_inplace<false,FrontValue,suffix_or_prefix_match_len>
                               (front,
                                d,
                                useq,
                                ulen,
                                vseq,
                                vlen,
                                seqnum0,
                                seqnum1);
    if ((vlen > ulen && vlen - ulen <= d) || (vlen <= ulen && ulen - vlen <= d))
    {
      /* Let vlen > ulen, then
            d+vlen >= ulen, as d is not negative. Moreover vlen - ulen <= d
            which implies d + vlen - ulen <= 2 * d.
         Len vlen <= ulen. Then from the above we know that ulen - vlen <= d
         which implies ulen <= d + vlen. Moreover, vlen - ulen < 0 and so
         d+vlen - ulen < d < 2 * d. So in both case d+vlen-ulen is in the
         range from 0 to 2*d and thus front[d+vlen-ulen] is defined. */
      if (front[d + vlen - ulen] >= ulen)
      {
        break;
      }
    }
    assert(d <= ulen + vlen);
    if constexpr (d_max_defined)
    {
      if (d == d_max)
      {
        free(front);
        return d_max+1;
      }
    }
  }
  free(front);
  return d;
}

template<typename char_type,class FrontValue,bool d_max_defined,
         LcsLcpLenType suffix_or_prefix_match_len>
static size_t fastedist_inplace(size_t d_max,
                                const char_type *useq,
                                size_t ulen,
                                const char_type *vseq,
                                size_t vlen,
                                size_t seqnum0,
                                size_t seqnum1)
{
  assert(ulen > 0 && vlen > 0);
  const size_t lcp_value0 = suffix_or_prefix_match_len(useq,0,vseq,0,
                                                       ulen,vlen,
                                                       seqnum0,seqnum1);

  if (ulen == vlen && lcp_value0 == ulen)
  {
    return 0;
  }
  FrontValue front0{lcp_value0};
  return fastedist_inplace_continue<char_type,FrontValue,0,d_max_defined,
                                    suffix_or_prefix_match_len>
                                   (&front0,
                                    d_max,
                                    useq,
                                    ulen,
                                    vseq,
                                    vlen,
                                    seqnum0,
                                    seqnum1);
}
#endif
