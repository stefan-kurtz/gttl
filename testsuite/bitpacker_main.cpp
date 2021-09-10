#include <cstddef>
#include <cstring>
#include <random>
#include <iostream>
#include <iomanip>
#include "utilities/runtime_class.hpp"
#include "utilities/mathsupport.hpp"
#include "utilities/unused.hpp"
#include "utilities/bitpacker.hpp"
#include "uint64_encoding.hpp"

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

#define COMPARE(V1,V2) if ((V1) != (V2))\
                       {\
                         std::cerr << #V1 << "=" << V1 << " != " \
                                   << V2 << "=" << #V2 << std::endl;\
                         exit(EXIT_FAILURE);\
                       }

#define RUN_TEST_CASES(COUNTER)\
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
            COMPARE(first_value,first_value_dec)\
            COMPARE(second_value,second_value_dec)\
            (COUNTER)++;\
          }\
        }

static void runner(bool direct,bool large,size_t num_values)
{
  size_t successes8 = 0, direct_successes8 = 0, successes9 = 0;
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
        if (large)
        {
          constexpr const int sizeof_unit = 9;
          GttlBitPacker<sizeof_unit,2> bp({first_bits,second_bits});
          RUN_TEST_CASES(successes9)
        }
      } else
      {
        if (!large)
        {
          if (direct)
          {
            Uint64Encoding<2> bp({first_bits,second_bits});
            for (size_t idx = 0; idx < num_values; idx++)
            {
              const uint64_t first_value = dis_first(rgen_first);
              for (uint64_t second_value = 0; second_value <= second_max;
                   second_value++)
              {
                const uint64_t code = bp.encode({first_value,second_value});
                const uint64_t first_value_dec = bp.decode_at<0>(code);
                const uint64_t second_value_dec = bp.decode_at<1>(code);
                COMPARE(first_value,first_value_dec)
                COMPARE(second_value,second_value_dec)
                direct_successes8++;
              }
            }
          } else
          {
            constexpr const int sizeof_unit = 8;
            GttlBitPacker<sizeof_unit,2> bp({first_bits,second_bits});
            RUN_TEST_CASES(successes8)
          }
        }
      }
    }
  }
  if (large)
  {
    std::cout << "# bitpack.successes9\t" << successes9 << std::endl;
  } else
  {
    if (direct)
    {
      std::cout << "# bitpack.direct_successes8\t" << direct_successes8
                << std::endl;
    } else
    {
      std::cout << "# bitpack.successes8\t" << successes8 << std::endl;
    }
  }
}

int main(int argc,char *argv[])
{
  long read_long;
  if (argc != 2 || sscanf(argv[1],"%ld",&read_long) != 1 || read_long <= 0)
  {
    std::cerr << "Usage: " << argv[0] << ": <num_of_value>" << std::endl;
    return EXIT_FAILURE;
  }
  const size_t num_values = static_cast<size_t>(read_long);
  RunTimeClass rt_small_direct;
  runner(true,false,num_values);
  rt_small_direct.show("8 bytes direct");
  RunTimeClass rt_small;
  runner(false,false,num_values);
  rt_small.show("8 bytes");
  RunTimeClass rt_large;
  runner(false,true,num_values);
  rt_small.show("9 bytes");
  return EXIT_SUCCESS;
}
