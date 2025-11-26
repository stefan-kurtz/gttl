#include <climits>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <random>
#include <iostream>
#include <format>
#include "utilities/bitpacker.hpp"
#include "utilities/runtime_class.hpp"
#include "utilities/mathsupport.hpp"
#include "utilities/bytes_unit.hpp"
#include "utilities/constexpr_for.hpp"
#include "uint64_encoding.hpp"

static void show_uint64_t_bytes([[maybe_unused]] uint64_t value)
{
#undef SHOWUINT64
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
                                   << V2 << "=" << #V2 << '\n';\
                         exit(EXIT_FAILURE);\
                       }

#define RUN_TEST_CASES\
        for (size_t idx = 0; idx < num_values; idx++)\
        {\
          const uint64_t first_value = dis_first(rgen_first);\
          for (uint64_t second_value = 0; second_value <= second_max; \
               second_value++)\
          {\
            show_uint64_t_bytes(first_value);\
            show_uint64_t_bytes(second_value);\
            const BytesUnit<sizeof_unit,2> \
                     bu(bp,{first_value,second_value});\
            const uint64_t first_value_dec = bu.template decode_at<0>(bp);\
            const uint64_t second_value_dec = bu.template decode_at<1>(bp);\
            COMPARE(first_value,first_value_dec)\
            COMPARE(second_value,second_value_dec)\
            successes++;\
          }\
        }

template<typename basetype,int overflow>
static void runner([[maybe_unused]] bool direct,size_t num_values)
{
  size_t successes = 0;
  std::mt19937_64 rgen_first;
  const int bits_basetype = sizeof(basetype) * CHAR_BIT;
  const int start_bits    = std::max(
                               9, static_cast<int>(bits_basetype / size_t(4)));
  for (int second_bits = 1; second_bits <= start_bits; second_bits++)
  {
    const uint64_t second_max = gttl_bits2maxvalue<uint64_t>(second_bits);

    for (int first_bits = 1; first_bits <= bits_basetype &&
                             first_bits + second_bits <= bits_basetype +
                                                         overflow * CHAR_BIT;
         first_bits++)
    {
      //std::cout << "first_bits\t" << first_bits
                //<< "\tsecond_bits\t" << second_bits << std::endl;
      const uint64_t first_max = gttl_bits2maxvalue<uint64_t>(first_bits);
      std::uniform_int_distribution<uint64_t> dis_first(0, first_max);
      if (first_bits + second_bits > bits_basetype)
      {
        if constexpr (overflow > 0)
        {
          static_assert(overflow <= 7);
          constexpr const int sizeof_unit = sizeof(basetype) + overflow;
          const GttlBitPacker<sizeof_unit, 2> bp({first_bits, second_bits});
          RUN_TEST_CASES
        }
      } else
      {
        if constexpr (overflow == 0)
        {
          if (direct)
          {
            if constexpr (sizeof(basetype) == sizeof(uint64_t))
            {
              const Uint64Encoding<2> bp({first_bits, second_bits});
              for (size_t idx = 0; idx < num_values; idx++)
              {
                const uint64_t first_value = dis_first(rgen_first);
                for (uint64_t second_value = 0; second_value <= second_max;
                     second_value++)
                {
                  const uint64_t code = bp.encode({first_value,second_value});
                  const uint64_t first_value_dec
                    = bp.template decode_at<0>(code);
                  const uint64_t second_value_dec
                    = bp.template decode_at<1>(code);
                  COMPARE(first_value,first_value_dec)
                  COMPARE(second_value,second_value_dec)
                  successes++;
                }
              }
            }
          } else
          {
            constexpr const int sizeof_unit = sizeof(basetype);
            const GttlBitPacker<sizeof_unit, 2> bp({first_bits, second_bits});
            RUN_TEST_CASES
          }
        }
      }
    }
  }
  std::cout << "# bitpacker.successes\t" << (sizeof(basetype) + overflow)
            << successes << '\n';
}

int main(int argc,char *argv[])
{
  long read_long;
  if (argc != 2 || sscanf(argv[1],"%ld",&read_long) != 1 || read_long <= 0)
  {
    std::cerr << "Usage: " << argv[0] << ": <num_of_value>\n";
    return EXIT_FAILURE;
  }
  const size_t num_values = static_cast<size_t>(read_long);

  for (int direct = 0; direct < 2; direct++)
  {
    RunTimeClass rt;
    runner<uint64_t,0>(direct ? true : false,num_values);
    rt.show(direct ? "8 bytes direct" : "8 bytes");
  }

  constexpr_for<0,7+1,1>([&](auto overflow)
  {
    RunTimeClass rt64;
    runner<uint64_t,overflow>(false,num_values);
    rt64.show(std::format("{} bytes", sizeof(uint64_t) + overflow));
    if constexpr (overflow <= 3)
    {
      RunTimeClass rt32;
      runner<uint32_t,overflow>(false,num_values);
      rt32.show(std::format("{} bytes", sizeof(uint32_t) + overflow));
    }
    if constexpr (overflow <= 1)
    {
      RunTimeClass rt16;
      runner<uint16_t,overflow>(false,num_values);
      rt16.show(std::format("{} bytes", sizeof(uint16_t) + overflow));
    }
  });
  return EXIT_SUCCESS;
}
