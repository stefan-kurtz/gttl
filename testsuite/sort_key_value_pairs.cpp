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
  union
  {
    double d;
    unsigned char ds[sizeof(double)];
  } key;
  size_t value;
  public:
  KeyValuePair(double _key,size_t _value)
  {
    this->key.d = _key;
    this->value = _value;
  }
  double double_get(void) const noexcept
  {
    return key.d;
  }
  void show(void) const noexcept
  {
    for (size_t idx = 0; idx < sizeof(double); idx++)
    {
      printf("%03d",static_cast<int>(key.ds[idx]));
      if (idx < sizeof(double) - 1)
      {
        printf(" ");
      }
    }
  }
  bool operator < (const KeyValuePair& other) const noexcept
  {
    return key.d < other.key.d;
  }
};

void show_key_values_pairs(const std::vector<KeyValuePair> &key_value_pairs)
{
  size_t idx = 0;
  for (auto &&kv : key_value_pairs)
  {
    printf("%03lu\t",idx);
    kv.show();
    printf("\n");
    idx++;
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
  size_t number_of_key_value_pairs = static_cast<size_t>(readlong);
  static_assert(sizeof(KeyValuePair) == 2 * sizeof(void *));
  std::vector<KeyValuePair> key_value_pairs{};
  key_value_pairs.reserve(number_of_key_value_pairs);
  UniformRandomDouble urd_gen(0,DBL_MAX);
  RunTimeClass rt_random_number_generation{};
  for (size_t idx = 0; idx < number_of_key_value_pairs; idx++)
  {
    double r = urd_gen.get();
    key_value_pairs.push_back(KeyValuePair(r,idx));
  }
  rt_random_number_generation.show("random number generation");
  std::cout << "# size of array with key value pairs (MB):\t"
            << (mega_bytes(key_value_pairs.size() * sizeof(KeyValuePair)))
            << std::endl;
  RunTimeClass rt_sorting{};
  if (use_radix_sort)
  {
    static constexpr const int sizeof_unit
      = static_cast<int>(sizeof(KeyValuePair));
    static constexpr const int num_sort_bits
      = static_cast<int>(CHAR_BIT * sizeof(double));
    const bool reversed_byte_order = is_big_endian() ? false : true;
    ska_large_lsb_small_radix_sort(sizeof_unit,
                                   num_sort_bits,
                                   reinterpret_cast<uint8_t *>
                                     (key_value_pairs.data()),
                                   key_value_pairs.size(),
                                   reversed_byte_order);
    StrFormat msg("sort %lu values with ska_large_lsb_small_radix_sort",
                  key_value_pairs.size());
    rt_sorting.show(msg.str());
  } else
  {
    std::sort(key_value_pairs.begin(),key_value_pairs.end());
    StrFormat msg("sort %lu values with std::sort",key_value_pairs.size());
    rt_sorting.show(msg.str());
  }
#ifdef SHOW_VALUES
  show_key_values_pairs(key_value_pairs);
#endif
  double previous = DBL_MIN;
  for (auto &&kv : key_value_pairs)
  {
    if (previous > kv.double_get())
    {
      std::cerr << argv[0] << ": previous=" << previous << " > "
                << kv.double_get() << "=next" << std::endl;
      return EXIT_FAILURE;
    }
    previous = kv.double_get();
  }
  return EXIT_SUCCESS;
}
