#ifndef IS_BIG_ENDIAN_HPP
#define IS_BIG_ENDIAN_HPP
#include <cstddef>
#include <cstdbool>

static inline bool is_big_endian(void)
{
  union
  {
     uint32_t integer;
     char bytes[4];
  } value = {0x01020304};
  return value.bytes[0] == static_cast<char>(1);
}

#endif
