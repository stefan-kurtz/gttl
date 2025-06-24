#ifndef REVERSE_COMPLEMENT_HPP
#define REVERSE_COMPLEMENT_HPP

#include <cstdint>
#include <cstddef>
#include "sequences/complement_uint8.hpp"

template<uint8_t (&complement_function)(uint8_t)>
static inline uint8_t *gttl_reverse_complement_encoded(const uint8_t *sequence,
                                                       size_t len)
{
  uint8_t *rc_seq = new uint8_t [len];

  for (size_t idx = 0; idx < len; idx++)
  {
    const uint8_t cc = sequence[idx];

    rc_seq[len-1-idx] = complement_function(cc);
  }
  return rc_seq;
}

static inline uint8_t *gttl_reverse_complement_encoded(const uint8_t *sequence,
                                                       size_t len)
{
  return gttl_reverse_complement_encoded<complement_uint8_wc_remains>
                                        (sequence,len);
}

template<uint8_t (&complement_function)(uint8_t)>
static inline void gttl_reverse_complement_encoded_inplace(uint8_t *sequence,
                                                           size_t len)
{
  uint8_t *fwd;
  uint8_t *bck;

  for (fwd = sequence, bck = sequence + len -1; fwd <= bck; fwd++, bck--)
  {
    const uint8_t tmp = *fwd;

    *fwd = complement_function(*bck);
    *bck = complement_function(tmp);
  }
}

static inline void gttl_reverse_complement_encoded_inplace(uint8_t *sequence,
                                                           size_t len)
{
  return gttl_reverse_complement_encoded_inplace<complement_uint8_wc_remains>
                                                (sequence,len);
}
#endif
