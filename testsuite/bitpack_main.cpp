#include <cstddef>
#include <cstring>
#include <random>
#include <iostream>
#include <iomanip>
#include "utilities/mathsupport.hpp"
#include "utilities/unused.hpp"
#include "utilities/bitpack.hpp"

static void show_uint64_t_bytes(GTTL_UNUSED uint64_t value)
{
#ifdef SHOWUINT64
  int shift = 56;
  for (size_t j = 0; j < sizeof(value); j++)
  {
    std::cout << std::right << std::setw(4)
              << static_cast<int>((value >> shift) & uint64_t(255))
              << (j == sizeof(value) - 1 ? "\n" : " ");
    shift -= CHAR_BIT;
  }
#endif
}

#define RUN_TEST_CASES\
        for (size_t idx = 0; idx < num_values; idx++)\
        {\
          uint8_t be_tmp[sizeof_unit];\
          const uint64_t first_value = dis_first(rgen_first);\
          for (uint64_t second_value = 0; second_value <= second_max; \
               second_value++)\
          {\
            show_uint64_t_bytes(first_value);\
            show_uint64_t_bytes(second_value);\
            const uint8_t *be = bp.encode({first_value,second_value});\
            memcpy(be_tmp,be,sizeof_unit);\
            const uint64_t first_value_dec = bp.decode_at<0>(be_tmp);\
            const uint64_t second_value_dec = bp.decode_at<1>(be_tmp);\
            if (first_value != first_value_dec)\
            {\
              std::cerr << "first_value=" << first_value << " != " \
                        << first_value_dec << "=first_value_dec" << std::endl;\
              exit(EXIT_FAILURE);\
            }\
            if (second_value != second_value_dec)\
            {\
              std::cerr << "second_value=" << second_value << " != " \
                        << second_value_dec << "=second_value_dec" \
                        << std::endl;\
              exit(EXIT_FAILURE);\
            }\
            successes++;\
          } \
        }

int main(int argc,char *argv[])
{
  long read_long;
  size_t successes = 0;
  if (argc != 2 || sscanf(argv[1],"%ld",&read_long) != 1 || read_long <= 0)
  {
    std::cerr << "Usage: " << argv[0] << ": <num_of_value>" << std::endl;
    return EXIT_FAILURE;
  }
  size_t num_values = static_cast<size_t>(read_long);
  std::mt19937_64 rgen_first;
  for (int second_bits = 1; second_bits <= 16; second_bits++)
  {
    const uint64_t second_max = GTTL_BITS2MAXVALUE(second_bits);

    for (int first_bits = 1; first_bits <= 64 &&
                             first_bits + second_bits <= 9 * CHAR_BIT;
         first_bits++)
    {
      const uint64_t first_max = GTTL_BITS2MAXVALUE(first_bits);
      std::uniform_int_distribution<uint64_t> dis_first(0, first_max);
      if (first_bits + second_bits > 64)
      {
        constexpr const int sizeof_unit = 9;
        GttlBitPack<sizeof_unit,2> bp({first_bits,second_bits});
        RUN_TEST_CASES
      } else
      {
        constexpr const int sizeof_unit = 8;
        GttlBitPack<sizeof_unit,2> bp({first_bits,second_bits});
        RUN_TEST_CASES
      }
    }
  }
  std::cout << "# bitpack.successes\t" << successes << std::endl;
  return EXIT_SUCCESS;
}
