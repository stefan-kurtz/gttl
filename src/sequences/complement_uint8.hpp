#ifndef COMPLEMENT_UINT8
#define COMPLEMENT_UINT8
#include <cstdint>

static inline uint8_t complement_uint8(uint8_t cc)
{
  return cc < uint8_t(4) ? (uint8_t(3) - cc) : cc;
}
#endif
