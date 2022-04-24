#include <cstdio>
#include <iostream>
#include <random>
#include <cfloat>
#include <tuple>
#include <vector>
#include <cstdio>
#include <iomanip>
#include "utilities/runtime_class.hpp"
#include "utilities/ska_lsb_radix_sort.hpp"
#include "utilities/is_big_endian.hpp"
#include "utilities/mathsupport.hpp"

class UniformRandomDouble
{
  std::random_device _rd{};
  std::mt19937 _gen{_rd()};
  std::uniform_real_distribution<double> _dist;

  void set(double low, double high)
  {
    std::uniform_real_distribution<double>::param_type param(low, high);
    _dist.param(param);
  }
  public:
  UniformRandomDouble(void)
  {
    set(DBL_MIN, DBL_MAX);
  }
  UniformRandomDouble(double low, double high)
  {
    set(low, high);
  }
  double get(void)
  {
    return _dist(_gen);
  }
};

class KeyValuePair
{
  union {
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
  double key_get(void) const noexcept
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
};

int main(int argc, char *argv[])
{
  long readlong;

  if (argc != 2 || sscanf(argv[1], "%ld", &readlong) != 1 || readlong <= 0)
  {
    std::cerr << "Usage: " << argv[0]
              << " <number of key/value (double/size_t>  pairs>" << std::endl;
    return EXIT_FAILURE;
  }
  size_t number_of_key_value_pairs = static_cast<size_t>(readlong);
  static_assert(sizeof(KeyValuePair) == 2 * sizeof(void *));
  std::vector<KeyValuePair> key_value_pairs{};
  key_value_pairs.reserve(number_of_key_value_pairs);
  UniformRandomDouble urd_gen{};
  RunTimeClass rt_random_number_generation{};
  for (size_t idx = 0; idx < number_of_key_value_pairs; idx++)
  {
    double r = urd_gen.get();
    key_value_pairs.push_back(KeyValuePair(r,idx));
  }
  rt_random_number_generation.show("random number generation");
  std::cout << "# size of array with key value pairs (MB):\t" <<
                  mega_bytes(key_value_pairs.size() * sizeof(KeyValuePair))
            << std::endl;
  RunTimeClass rt_sorting{};
  static constexpr const int sizeof_unit
    = static_cast<int>(sizeof(KeyValuePair));
  static constexpr const int num_sort_bits
    = static_cast<int>(CHAR_BIT * sizeof(double));
  const bool reversed_byte_order = is_big_endian() ? false : true;
  using BucketCounttype = size_t;
  Buckets<BucketCounttype> *buckets
    = ska_lsb_radix_sort<BucketCounttype>(sizeof_unit,
                                          num_sort_bits,
                                          reinterpret_cast<uint8_t *>
                                            (key_value_pairs.data()),
                                          key_value_pairs.size(),
                                          reversed_byte_order);
  size_t idx = 0;
  for (auto it = buckets->begin(); it != buckets->end(); ++it)
  {
    std::cout << "# bucket\t" << idx << "\t"
              << (std::get<1>(*it) - std::get<0>(*it))
              << std::endl;
    idx++;
  }
  std::cout << "# maximum width of buckets\t"
            << buckets->maximum_width_get() << std::endl;
  delete buckets;
  rt_sorting.show("sorting");
  double previous = DBL_MIN;
  for (auto &&kv : key_value_pairs)
  {
    assert(previous <= kv.key_get());
#define SHOW
#ifdef SHOW
    kv.show();
    std::cout << std::endl;
#endif
    previous = kv.key_get();
  }
  return EXIT_SUCCESS;
}
