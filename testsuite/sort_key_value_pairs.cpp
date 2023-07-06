#include <cstdio>
#include <tuple>
#include <vector>
#include <cstdio>
#include <cstring>
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

static std::string int2byte(int i)
{
  assert(i <= UINT8_MAX and i >= 0);
  const size_t leading_zeros = i < 10 ? 2 : (i < 100 ? 1 : 0);
  return std::string(leading_zeros,'0') + std::to_string(i);
}

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
    const uint8_t *ds = reinterpret_cast<const uint8_t *>(&key);
    for (size_t idx = 0; idx < sizeof(double); idx++)
    {
      s += int2byte(static_cast<int>(ds[idx]));
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

class Key2ValuePair
{
  static constexpr const size_t sizeof_key = 2 * sizeof(double);
  uint8_t keys_as_bytes[sizeof_key];
  size_t value;
  public:
  Key2ValuePair(void)
    : value(0)
  {
    memset(keys_as_bytes,0,sizeof_key);
  }
  Key2ValuePair(double key0,double key1,size_t _value)
    : value(_value)
  {
    const uint8_t *ds = reinterpret_cast<const uint8_t *>(&key0);
    for (size_t idx = 0; idx < sizeof(double); idx++)
    {
      keys_as_bytes[idx] = ds[idx];
    }
    ds = reinterpret_cast<const uint8_t *>(&key1);
    for (size_t idx = 0; idx < sizeof(double); idx++)
    {
      keys_as_bytes[sizeof(double) + idx] = ds[idx];
    }
  }
  std::string to_string(void) const noexcept
  {
    std::string s{};
    for (size_t idx = 0; idx < sizeof_key; idx++)
    {
      s += int2byte(static_cast<int>(keys_as_bytes[idx]));
      if (idx < sizeof_key - 1)
      {
        s += " ";
      }
    }
    s += '\t' + std::to_string(value);
    return s;
  }
  bool operator < (const Key2ValuePair& other) const noexcept
  {
    return memcmp(keys_as_bytes,other.keys_as_bytes,sizeof_key) < 0;
  }
};

template<class T>
static void sort_values(unsigned int seed,
                        const char *progname,
                        bool show,
                        bool use_radix_sort,
                        size_t number_of_values,
                        double max_random,
                        int num_sort_bits)
{
  std::vector<T> values{};
  values.reserve(number_of_values);
  UniformRandomDouble urd_gen(0,max_random,seed);
  RunTimeClass rt_random_number_generation{};
  for (size_t idx = 0; idx < number_of_values; idx++)
  {
    const double r = urd_gen.get();
    if constexpr (std::is_same_v<T, KeyValuePair>)
    {
      values.push_back(T(r,idx));
    } else
    {
      if constexpr (std::is_same_v<T, Key2ValuePair>)
      {
        const double s = urd_gen.get();
        values.push_back(T(r,s,idx));
      } else
      {
        values.push_back(static_cast<T>(r));
      }
    }
  }
  rt_random_number_generation.show("random number generation");
  std::cout << "# size of array (MB):\t"
            << (mega_bytes(values.size() * sizeof(T)))
            << std::endl;
  RunTimeClass rt_sorting{};
  const char *tag;
  if constexpr (std::is_same_v<T, KeyValuePair>)
  {
    tag = "key value pairs";
  } else
  {
    if constexpr (std::is_same_v<T, Key2ValuePair>)
    {
      tag = "key0+key1 value pairs";
    } else
    {
      tag = "integers";
    }
  }
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
      if constexpr (std::is_same_v<T, Key2ValuePair>)
      {
        static constexpr const int sizeof_unit
          = static_cast<int>(sizeof(T));
        ska_large_lsb_small_radix_sort(sizeof_unit,
                                       num_sort_bits,
                                       reinterpret_cast<uint8_t *>
                                                       (values.data()),
                                       values.size(),
                                       false);
      } else
      {
        static_assert(sizeof(T) == sizeof(void *));
        ska_large_lsb_small_radix_sort(num_sort_bits,
                                       reinterpret_cast<uint64_t *>
                                                       (values.data()),
                                       values.size());
      }
    }
    StrFormat msg("sort %lu %s with ska_large_lsb_small_radix_sort",
                  values.size(),tag);
    rt_sorting.show(msg.str());
  } else
  {
    std::sort(values.begin(),values.end());
    StrFormat msg("sort %lu %s with std::sort",values.size(),tag);
    rt_sorting.show(msg.str());
  }
  if (show)
  {
    for (auto v : values)
    {
      if constexpr (std::is_same_v<T, KeyValuePair> ||
                    std::is_same_v<T, Key2ValuePair>)
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
      if constexpr (std::is_same_v<T, KeyValuePair> ||
                    std::is_same_v<T, Key2ValuePair>)
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
      (strcmp(argv[1],"p") != 0 && strcmp(argv[1],"i") != 0 &&
       strcmp(argv[1],"t") != 0) ||
      (strcmp(argv[2],"rad") != 0 && strcmp(argv[2],"std") != 0) ||
       sscanf(argv[3], "%ld", &readlong) != 1 || readlong <= 0)
  {
    std::cerr << "Usage: " << argv[0] << " p|i rad|std <number of values>"
              << std::endl << "i     create array of integers"
              << std::endl << "p     create array of pairs"
              << std::endl << "t     create array of triples (fst 2 are keys)"
              << std::endl << "rad   use radix sort"
              << std::endl << "std   use std::sort"
              << std::endl;
    return EXIT_FAILURE;
  }
  const bool use_radix_sort = strcmp(argv[2],"rad") == 0 ? true : false;
  const size_t number_of_values = static_cast<size_t>(readlong);
  const bool show = false;
  unsigned int seed = 0; /* use some fixed value > 0 for a fixed seed */
  if (strcmp(argv[1],"p") == 0)
  {
    static constexpr const int num_sort_bits
      = static_cast<int>(CHAR_BIT * sizeof(double));
    sort_values<KeyValuePair>(seed,argv[0],show,use_radix_sort,number_of_values,
                              DBL_MAX,num_sort_bits);
  } else
  {
    if (strcmp(argv[1],"t") == 0)
    {
      static constexpr const int num_sort_bits
        = static_cast<int>(CHAR_BIT * 2 * sizeof(uint64_t));
      sort_values<Key2ValuePair>(seed,argv[0],show,use_radix_sort,
                                 number_of_values,
                                 DBL_MAX,num_sort_bits);
    } else
    {
      //gttl_required_bits<size_t>(number_of_values);
      const int num_sort_bits = 64;
      sort_values<size_t>(seed,argv[0],show,use_radix_sort,number_of_values,
                          static_cast<double>(number_of_values),
                          num_sort_bits);
    }
  }
  return EXIT_SUCCESS;
}
