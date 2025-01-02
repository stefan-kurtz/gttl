#ifndef REVERSE_COMPLEMENT_HPP
#define REVERSE_COMPLEMENT_HPP

#include <cstdint>
#include <cstddef>
#include "sequences/complement_uint8.hpp"

static inline uint8_t *gttl_reverse_complement_encoded(const uint8_t *sequence,
                                                       size_t len)
{
  uint8_t *rc_seq = new uint8_t [len];

  for (size_t idx = 0; idx < len; idx++)
  {
    const uint8_t cc = sequence[idx];

    rc_seq[len-1-idx] = complement_uint8_wc_remains(cc);
  }
  return rc_seq;
}

static inline void gttl_reverse_complement_encoded_inplace(uint8_t *sequence,
                                                           size_t len)
{
  uint8_t *fwd, *bck;

  for (fwd = sequence, bck = sequence + len -1; fwd <= bck; fwd++, bck--)
  {
    const uint8_t tmp = *fwd;

    *fwd = complement_uint8_wc_remains(*bck);
    *bck = complement_uint8_wc_remains(tmp);
  }
}
#endif
