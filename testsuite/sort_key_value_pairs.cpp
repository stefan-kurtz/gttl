#include <cstdio>
#include <tuple>
#include <vector>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <type_traits>
#include <algorithm>
#include "utilities/runtime_class.hpp"
#include "utilities/mathsupport.hpp"
#include "utilities/str_format.hpp"
#include "utilities/is_big_endian.hpp"
#include "utilities/uniform_random_double.hpp"
#include "utilities/ska_lsb_radix_sort.hpp"

class KeyValuePair
{
  double key;
  size_t value;
  public:
  KeyValuePair(void)
    : key(0)
    , value(0)
  {}
  KeyValuePair(double _key,size_t _value)
    : key(_key)
    , value(_value)
  {}
  std::string to_string(void) const noexcept
  {
    std::string s{};
    const unsigned char *ds = reinterpret_cast<const unsigned char *>(&key);
    for (size_t idx = 0; idx < sizeof(double); idx++)
    {
      char sbuf[3+1];
      sprintf(sbuf,"%03d",static_cast<int>(ds[idx]));
      s += std::string(sbuf);
      if (idx < sizeof(double) - 1)
      {
        s += " ";
      }
    }
    s += '\t' + std::to_string(value);
    return s;
  }
  bool operator < (const KeyValuePair& other) const noexcept
  {
    return key < other.key;
  }
};

template<class T>
static void sort_values(unsigned int seed,
                        const char *progname,bool show,bool use_radix_sort,
                        size_t number_of_values,double max_random,
                        int num_sort_bits)
{
  std::vector<T> values{};
  values.reserve(number_of_values);
  UniformRandomDouble urd_gen(0,max_random,seed);
  RunTimeClass rt_random_number_generation{};
  for (size_t idx = 0; idx < number_of_values; idx++)
  {
    double r = urd_gen.get();
    if constexpr (std::is_same_v<T, KeyValuePair>)
    {
      values.push_back(T(r,idx));
    } else
    {
      values.push_back(static_cast<T>(r));
    }
  }
  rt_random_number_generation.show("random number generation");
  std::cout << "# size of array (MB):\t"
            << (mega_bytes(values.size() * sizeof(T)))
            << std::endl;
  RunTimeClass rt_sorting{};
  if (use_radix_sort)
  {
    if constexpr (std::is_same_v<T, KeyValuePair>)
    {
      static constexpr const int sizeof_unit
        = static_cast<int>(sizeof(T));
      const bool reversed_byte_order = is_big_endian() ? false : true;
      ska_large_lsb_small_radix_sort(sizeof_unit,
                                     num_sort_bits,
                                     reinterpret_cast<uint8_t *>(values.data()),
                                     values.size(),
                                     reversed_byte_order);
    } else
    {
      static_assert(sizeof(T) == sizeof(void *));
      ska_large_lsb_small_radix_sort(num_sort_bits,
                                     reinterpret_cast<uint64_t *>
                                                     (values.data()),
                                     values.size());
    }
    StrFormat msg("sort %lu values with ska_large_lsb_small_radix_sort",
                  values.size());
    rt_sorting.show(msg.str());
  } else
  {
    std::sort(values.begin(),values.end());
    StrFormat msg("sort %lu values with std::sort",values.size());
    rt_sorting.show(msg.str());
  }
  if (show)
  {
    for (auto v : values)
    {
      if constexpr (std::is_same_v<T, KeyValuePair>)
      {
        std::cout << v.to_string() << std::endl;
      } else
      {
        std::cout << v << std::endl;
      }
    }
  }
  T previous;
  bool previous_defined = false;
  for (auto v : values)
  {
    if (previous_defined && v < previous)
    {
      if constexpr (sizeof(T) == 2 * sizeof(void *))
      {
        std::cerr << progname << ": previous=" << previous.to_string()
                  << " > " << v.to_string() << "=next" << std::endl;
      } else
      {
        std::cerr << progname << ": previous=" << previous
                  << " > " << v << "=next" << std::endl;
      }
      exit(EXIT_FAILURE);
    }
    previous = v;
    previous_defined = true;
  }
}

int main(int argc, char *argv[])
{
  static_assert(sizeof(KeyValuePair) == 2 * sizeof(void *));
  long readlong;

  if (argc != 4 ||
      (strcmp(argv[1],"p") != 0 && strcmp(argv[1],"i") != 0) ||
      (strcmp(argv[2],"rad") != 0 && strcmp(argv[2],"std") != 0) ||
       sscanf(argv[3], "%ld", &readlong) != 1 || readlong <= 0)
  {
    std::cerr << "Usage: " << argv[0] << " p|i rad|std <number of values>"
              << std::endl << "p     create array of pairs"
              << std::endl << "i     create array of integers"
              << std::endl << "rad   use radix sort"
              << std::endl << "std   use std::sort"
              << std::endl;
    return EXIT_FAILURE;
  }
  const bool use_pairs = strcmp(argv[1],"p") == 0 ? true : false;
  const bool use_radix_sort = strcmp(argv[2],"rad") == 0 ? true : false;
  const size_t number_of_values = static_cast<size_t>(readlong);
  const bool show = false;
  unsigned int seed = 0; /* use some fixed value > 0 for a fixed seed */
  if (use_pairs)
  {
    static constexpr const int num_sort_bits
      = static_cast<int>(CHAR_BIT * sizeof(double));
    sort_values<KeyValuePair>(seed,argv[0],show,use_radix_sort,number_of_values,
                              DBL_MAX,num_sort_bits);
  } else
  {
    const int num_sort_bits = 64;//gttl_required_bits<size_t>(number_of_values);
    sort_values<size_t>(seed,argv[0],show,use_radix_sort,number_of_values,
                        static_cast<double>(number_of_values),
                        num_sort_bits);
  }
  return EXIT_SUCCESS;
}
