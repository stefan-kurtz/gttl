#include <cstdio>
#include <tuple>
#include <vector>
#include <cstdio>
#include <iomanip>
#include <iostream>
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
  double key_get(void) const noexcept
  {
    return key;
  }
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
    return s;
  }
  bool operator < (const KeyValuePair& other) const noexcept
  {
    return key < other.key;
  }
};

template<class T>
static void sort_values(const char *progname,bool use_radix_sort,
                        size_t number_of_values,
                        int num_sort_bits)
{
  std::vector<T> values{};
  values.reserve(number_of_values);
  UniformRandomDouble urd_gen(0,DBL_MAX);
  RunTimeClass rt_random_number_generation{};
  for (size_t idx = 0; idx < number_of_values; idx++)
  {
    double r = urd_gen.get();
    values.push_back(T(r,idx));
  }
  rt_random_number_generation.show("random number generation");
  std::cout << "# size of array (MB):\t"
            << (mega_bytes(values.size() * sizeof(T)))
            << std::endl;
  RunTimeClass rt_sorting{};
  if (use_radix_sort)
  {
    static constexpr const int sizeof_unit
      = static_cast<int>(sizeof(T));
    const bool reversed_byte_order = is_big_endian() ? false : true;
    ska_large_lsb_small_radix_sort(sizeof_unit,
                                   num_sort_bits,
                                   reinterpret_cast<uint8_t *>(values.data()),
                                   values.size(),
                                   reversed_byte_order);
    StrFormat msg("sort %lu values with ska_large_lsb_small_radix_sort",
                  values.size());
    rt_sorting.show(msg.str());
  } else
  {
    std::sort(values.begin(),values.end());
    StrFormat msg("sort %lu values with std::sort",values.size());
    rt_sorting.show(msg.str());
  }
#ifdef SHOW_VALUES
  for (auto v : values)
  {
    std::cout << v.to_string() << std::endl;
  }
#endif
  T previous;
  bool previous_defined = false;
  for (auto v : values)
  {
    if (previous_defined && v < previous)
    {
      std::cerr << progname << ": previous=" << previous.to_string() << " > "
                << v.to_string() << "=next" << std::endl;
      exit(EXIT_FAILURE);
    }
    previous = v;
    previous_defined = true;
  }
}

int main(int argc, char *argv[])
{
  long readlong;

  if (argc != 3 || (strcmp(argv[1],"std") != 0 && strcmp(argv[1],"rad") != 0) ||
      sscanf(argv[2], "%ld", &readlong) != 1 || readlong <= 0)
  {
    std::cerr << "Usage: " << argv[0]
              << "rad|std <number of key/value (double/size_t) pairs>"
              << std::endl << "rad triggers the use of radix sort"
              << std::endl << "std triggers the use of std::sort"
              << std::endl;
    return EXIT_FAILURE;
  }
  const bool use_radix_sort = strcmp(argv[1],"rad") == 0 ? true : false;
  const size_t number_of_values = static_cast<size_t>(readlong);
  static_assert(sizeof(KeyValuePair) == 2 * sizeof(void *));
  static constexpr const int num_sort_bits
    = static_cast<int>(CHAR_BIT * sizeof(double));
  sort_values<KeyValuePair>(argv[0],use_radix_sort,number_of_values,
                            num_sort_bits);
  return EXIT_SUCCESS;
}
