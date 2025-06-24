#ifndef IS_BIG_ENDIAN_HPP
#define IS_BIG_ENDIAN_HPP
#include <bit>

#if __cplusplus > 201703L
consteval bool is_big_endian(void)
{
  if constexpr (std::endian::native == std::endian::big)
  {
    return true;
  } else
  {
    static_assert(std::endian::native == std::endian::little);
    return false;
  }
}
#else
/* clang++ does not allow this function to be a constexpr. So we omit it.  */
bool is_big_endian(void)
{
  union
  {
     uint32_t integer;
     char bytes[4];
  } value{0x01020304};
  return value.bytes[0] == static_cast<char>(1);
}
#endif
#endif
