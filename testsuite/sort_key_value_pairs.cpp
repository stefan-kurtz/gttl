#include <cstdio>
#include <tuple>
#include <vector>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <type_traits>
#include <algorithm>
#include <iterator>
#include "utilities/cxxopts.hpp"
#include "utilities/runtime_class.hpp"
#include "utilities/mathsupport.hpp"
#include "utilities/str_format.hpp"
#include "utilities/is_big_endian.hpp"
#include "utilities/uniform_random_double.hpp"
#include "utilities/ska_lsb_radix_sort.hpp"
#include "utilities/merge_sort.hpp"
#include "utilities/runtime_class.hpp"

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << std::endl;
}

class SortKeyValuePairsOptions
{
 private:
  size_t number_of_values;
  size_t num_threads;
  bool help_option;
  char data_type_option;
  int sort_mode;
  std::string sort_mode_option;

 public:
  SortKeyValuePairsOptions()
    : number_of_values(0)
    , num_threads(size_t(1))
    , data_type_option('i')
    , sort_mode(0)
    , sort_mode_option("lsb-radix")
  {}

  void parse(int argc, char **argv)
  {
    cxxopts::Options options(argv[0],"run sorting methods for different "
                                     "kinds of data");
    options.set_width(80);
    options.custom_help(std::string("[options] number_of_values"));
    options.set_tab_expansion();
    options.add_options()
       ("d,data_type",
        "specify the data type to be sorted: i means integers, "
        "p means pairs of integers (where the first integer is considered "
        "the key), t means triples of integers (where the first two integers "
        "are considered the key)",
        cxxopts::value<char>(data_type_option)->default_value("i"))
       ("m,sort_mode",
        "specify the method to use for sorting, "
        "possible are: lsb-radix, mergesort or stdsort",
        cxxopts::value<std::string>(sort_mode_option)
                                   ->default_value("lsb-radix"))
       ("t,num_threads",
        "specify the number of threads",
        cxxopts::value<size_t>(num_threads)->default_value("1"))
       ("h,help", "print usage");
    try
    {
      auto result = options.parse(argc, argv);
      if (result.count("help") > 0)
      {
        help_option = true;
        usage(options);
      }
      const std::vector<std::string>& unmatched_args = result.unmatched();
      if (unmatched_args.size() < 1)
      {
        throw std::invalid_argument("missing positional number_of_values "
                                    "argument");
      } else
      {
        if (unmatched_args.size() > 1)
        {
          throw std::invalid_argument("superfluous positional argument");
        } else
        {
          long read_long;
          if (std::sscanf(unmatched_args[0].c_str(),"%ld",&read_long) != 1 or
              read_long < 0)
          {
            throw std::invalid_argument("positional argument must be positive "
                                        "integer specifying the number of "
                                        "values to sort");
          }
          number_of_values = static_cast<size_t>(read_long);
        }
      }
      if (data_type_option != 'i' and data_type_option != 'p' and
          data_type_option != 't')
      {
        throw std::invalid_argument("argument to option -d/--data_type must be "
                                    "i, p, or t");
      }
      if (sort_mode_option == std::string("lsb-radix"))
      {
        sort_mode = 0;
        if (data_type_option != 'i' && num_threads != size_t(1))
        {
          throw std::invalid_argument("option -t/--num_threads can only be "
                                      "used for data_type i");
        }
      } else
      {
        if (num_threads != size_t(1))
        {
          throw std::invalid_argument("option -t/--num_threads can only be "
                                      "used for lsb-radix");
        }
        if (sort_mode_option == std::string("mergesort"))
        {
          sort_mode = 1;
        } else
        {
          if (sort_mode_option == std::string("stdsort"))
          {
            sort_mode = 2;
          } else
          {
            throw std::invalid_argument("argument to option -m/--sort_mode "
                                        "must be lsb-radix or mergesort or "
                                        "stdsort");
          }
        }
      }
    }
    catch (const cxxopts::OptionException &e)
    {
      usage(options);
      throw std::invalid_argument(e.what());
    }
  }
  size_t number_of_values_get(void) const noexcept
  {
    return number_of_values;
  }
  size_t num_threads_get(void) const noexcept
  {
    return num_threads;
  }
  char data_type_option_get(void) const noexcept
  {
    return data_type_option;
  }
  int sort_mode_get(void) const noexcept
  {
    return sort_mode;
  }
  bool help_option_is_set(void) const noexcept
  {
    return help_option;
  }
};

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

template<class It>
static inline void show_values(It begin, It end)
{
  for (auto it = begin; it != end; ++it)
  {
    auto v = *it;
    using T = decltype(v);
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

template<class T>
static void sort_values(unsigned int seed,
                        const char *progname,
                        bool show,
                        int sort_mode,
                        size_t number_of_values,
                        double max_random,
                        int num_sort_bits,
                        size_t num_threads)
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
  if (sort_mode == 0)
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
        static_assert(sizeof(T) == sizeof(uint64_t));
        ska_large_lsb_small_radix_sort(num_sort_bits,
                                       reinterpret_cast<uint64_t *>
                                                       (values.data()),
                                       values.size(),
                                       num_threads);
      }
    }
    StrFormat msg("sort %zu %s with ska_large_lsb_small_radix_sort and %zu "
                  "thread%s",
                  values.size(),tag,num_threads,num_threads == 1 ? "" : "s");
    rt_sorting.show(msg.str());
  } else
  {
    if (sort_mode == 1)
    {
      const unsigned int n_threads = 1;
      merge_sort<decltype(values.begin())>(values.begin(),values.end(),
                                           n_threads);
      StrFormat msg("sort %zu %s with mergesort",values.size(),tag);
      rt_sorting.show(msg.str());
    } else
    {
      std::sort(values.begin(),values.end());
      StrFormat msg("sort %zu %s with std::sort",values.size(),tag);
      rt_sorting.show(msg.str());
    }
  }
  if (show)
  {
    show_values<decltype(values.begin())>(values.begin(),values.end());
  }
  if (values.size() > 1)
  {
    T previous = values[0];
    for (size_t idx = 1; idx < values.size(); idx++)
    {
      T v = values[idx];
      if (v < previous)
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
    }
  }
}

int main(int argc, char *argv[])
{
  static_assert(sizeof(KeyValuePair) == 2 * sizeof(void *));

  SortKeyValuePairsOptions options;
  try
  {
    options.parse(argc,argv);
  }
  catch (std::invalid_argument &e) /* check_err.py */
  {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  if (options.help_option_is_set())
  {
    return EXIT_SUCCESS;
  }
  const bool show = false;
  const int sort_mode = options.sort_mode_get();
  unsigned int seed = 0; /* use some fixed value > 0 for a fixed seed */
  if (options.data_type_option_get() == 'p')
  {
    static constexpr const int num_sort_bits
      = static_cast<int>(CHAR_BIT * sizeof(double));
    sort_values<KeyValuePair>(seed,argv[0],show,sort_mode,
                              options.number_of_values_get(),
                              DBL_MAX,num_sort_bits,
                              size_t(1));
  } else
  {
    if (options.data_type_option_get() == 't')
    {
      static constexpr const int num_sort_bits
        = static_cast<int>(CHAR_BIT * 2 * sizeof(uint64_t));
      sort_values<Key2ValuePair>(seed,argv[0],show,sort_mode,
                                 options.number_of_values_get(),
                                 DBL_MAX,num_sort_bits,
                                 size_t(1));
    } else
    {
      assert(options.data_type_option_get() == 'i');
      constexpr const int num_sort_bits = 64;
      const double max_random
        = static_cast<double>(options.number_of_values_get());
      sort_values<uint64_t>(seed,argv[0],show,sort_mode,
                            options.number_of_values_get(),
                            max_random,
                            num_sort_bits,
                            options.num_threads_get());
    }
  }
  return EXIT_SUCCESS;
}
