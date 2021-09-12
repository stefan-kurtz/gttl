#include <climits>
#include <iostream>
#include "utilities/is_big_endian.hpp"

/*#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
  reinterpret cast does not work
#endif
*/

static void byte_order(void)
{
  uint8_t bytes[8];
  for (int idx = 0; idx < 8; idx++)
  {
    bytes[idx] = static_cast<uint8_t>(idx);
  }
  uint64_t integer = *(reinterpret_cast<uint64_t *>(bytes));
  int shift = 56;
  std::cout << "# byte index (left to right)\tinteger order" << std::endl;
  for (int idx = 0; idx < 8; idx++)
  {
    std::cout << idx << "\t"
              << static_cast<int>((integer >> shift) & uint64_t(255))
              << std::endl;
    shift -= CHAR_BIT;
  }
}

int main(void)
{
  std::cout << "is_big_endian\t" << (is_big_endian() ? "true" : "false")
            << std::endl;
  byte_order();
}
