#ifndef PREFIX_SUFFIX_EXTENSION_HPP
#define PREFIX_SUFFIX_EXTENSION_HPP
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>
#include "sequences/lcs_lcp_len_type.hpp"
#include "sequences/outsenseedist_inplace.hpp"

/* run the iteration with increasing value of d and return the
   edist distance, if the last row and last column was reached.
*/

template<bool track_eop,
         class Tracker,
         class FrontValue,
         LcsLcpLenType suffix_or_prefix_match_len>
static void prefix_or_suffix_extension_generic(Tracker *tracker,
                                               const char *useq,
                                               size_t ulen,
                                               const char *vseq,
                                               size_t vlen,
                                               size_t useqnum,
                                               size_t vseqnum)
{
  if (ulen == 0 || vlen == 0)
  {
    return;
  }
  const size_t initial_lcp
    = suffix_or_prefix_match_len(useq,0,vseq,0,ulen,vlen,useqnum,vseqnum);
  assert(initial_lcp <= ulen && initial_lcp <= vlen);
  std::vector<FrontValue> front{};
  size_t upper_bound_d = 32;
  front.reserve(2 * upper_bound_d + 1);
  front.push_back(FrontValue(initial_lcp));
  if (tracker->evaluate(0, 0, 0, front.data() + 0) == 0)
  {
    return;
  }
  for (size_t d = 1; /* Nothing */; d++)
  {
    if (d > upper_bound_d)
    {
      upper_bound_d += 32;
      front.reserve(2 * upper_bound_d + 1);
    }
    front.push_back(FrontValue());
    front.push_back(FrontValue());
    outsense_next_front_inplace<track_eop,FrontValue,suffix_or_prefix_match_len>
                               (front.data(),
                                d,
                                useq,
                                ulen,
                                vseq,
                                vlen,
                                useqnum,
                                vseqnum);
    if (tracker->evaluate(d, -static_cast<int32_t>(d), static_cast<int32_t>(d),
                          front.data() + d) == 0)
    {
      break;
    }
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
  }
}
#endif
