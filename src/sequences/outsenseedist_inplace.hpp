#ifndef OUTSENSEEDIST_INPLACE_HPP
#define OUTSENSEEDIST_INPLACE_HPP
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include "sequences/lcs_lcp_len_type.hpp"

template<class FrontValue,LcsLcpLenType suffix_or_prefix_match_len>
static inline void outsense_next_front_inplace(FrontValue *front,
                                               size_t d,
                                               const char *useq,
                                               size_t ulen,
                                               const char *vseq,
                                               size_t vlen,
                                               size_t seqnum0,
                                               size_t seqnum1)
{
  assert(d > 0 && front != nullptr);
  FrontValue insertion_value = front[0]; /* copy constructor */
  FrontValue lvalue = insertion_value; /* copy constructor */
  lvalue += 1; /* add_error, do not replace by ++ */
  /* front(-d,d) */
  front[0] = lvalue;
  if (lvalue < ulen && lvalue < vlen + d)
  {
    assert(lvalue + 0 >= d);
    const size_t l = suffix_or_prefix_match_len(useq,lvalue + 0,
                                                vseq,(lvalue + 0) - d,
                                                ulen,vlen,
                                                seqnum0,seqnum1);
    front[0] += l; /* add match of length l */
  }
  FrontValue replacement_value{};
  if (d > size_t(1))
  {
    replacement_value = front[1];
    if (lvalue <= replacement_value)
    {
      lvalue = replacement_value;
      lvalue += 1; /* do not replace by ++ */
    }
  }
  front[1] = lvalue;
  if (lvalue < ulen && lvalue + 1 < vlen + d)
  {
    assert(lvalue + 1 >= d);
    const size_t l = suffix_or_prefix_match_len(useq,lvalue + 0,
                                                vseq,(lvalue + 1) - d,
                                                ulen,vlen,
                                                seqnum0,seqnum1);
    front[1] += l;
  }
  const size_t frontmaxidx = 2 * d;
  for (size_t idx=size_t(2); idx <= frontmaxidx; idx++)
  {
    lvalue = insertion_value; /* insertion */
    lvalue += 0; /* do not delete */
    if (idx <= frontmaxidx - 1 && lvalue <= replacement_value)
    {
      lvalue = replacement_value; /* replacement */
      lvalue += 1;
    }
    if (idx <= frontmaxidx - 2 && lvalue <= front[idx])
    {
      lvalue = front[idx]; /* deletion */
      lvalue += 1;
    }
    insertion_value = replacement_value;
    if (idx < frontmaxidx)
    {
      replacement_value = front[idx];
    }
    front[idx] = lvalue;
    if (lvalue < ulen && lvalue + idx < vlen + d)
    {
      assert(lvalue + idx >= d);
      const size_t l = suffix_or_prefix_match_len(useq,lvalue + 0,
                                                  vseq,(lvalue + idx) - d,
                                                  ulen,vlen,
                                                  seqnum0,seqnum1);
      front[idx] += l;
    }
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
  size_t d, allocated = 0;

  if constexpr (d_max_defined)
  {
    assert(d_max > 0 && static_cast<size_t>(previous_d) < d_max);
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
    outsense_next_front_inplace<FrontValue,suffix_or_prefix_match_len>
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
