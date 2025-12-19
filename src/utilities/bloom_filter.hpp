/*
** Developed by Henning Lindemann
 */

#ifndef BLOOM_FILTER_HPP
#define BLOOM_FILTER_HPP

#include "utilities/bloom_filter_hash_function.hpp"
#include "utilities/bloom_filter_u64.hpp"
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <numbers>
#include <vector>

template <bool thread_safe = true>
class BloomFilter
{
  private:
  double ln_error; /* only used in one case */
  size_t num_hash_functions;
  size_t num_bits;
  std::vector<BloomFilterU64<thread_safe>> data_vec;
  bool set_bit(uint64_t idx)
  {
    return data_vec[idx / 64].set_bit(idx % 64);
  }
  [[nodiscard]] size_t num_bits_3args(double error,
                                      double d_number_of_elements,
                                      double d_num_hash_functions) const
  {
    return static_cast<size_t>
           ((-1.0 * d_number_of_elements * d_num_hash_functions) /
            std::log(1.0 - std::pow(error, 1.0/d_num_hash_functions)));
  }

  public:
  [[nodiscard]] bool get_bit(uint64_t index) const
  {
    return data_vec[index / 64].get_bit(index % 64);
  }

  BloomFilter(size_t _num_bits, size_t _num_hash_functions)
    : ln_error(0.0)
    , num_hash_functions(_num_hash_functions)
    , num_bits(_num_bits)
    , data_vec(std::vector<BloomFilterU64<thread_safe>>((num_bits + 63) / 64))
  { }

  BloomFilter(double error, size_t number_of_elements)
    // the numbers are calculated according to the recommendations
    // in https://en.wikipedia.org/wiki/Bloom_filter
    : ln_error(std::log(error))
    , num_hash_functions(static_cast<size_t>(ln_error * -std::numbers::log2e))
    , num_bits(static_cast<size_t>(static_cast<double>(number_of_elements) *
                                   ln_error * -2.081368))
    , data_vec(std::vector<BloomFilterU64<thread_safe>>((num_bits + 63) / 64))
  { }

  BloomFilter(double error, size_t number_of_elements,
              size_t _num_hash_functions)
    : ln_error(0.0)
    , num_hash_functions(_num_hash_functions)
    , num_bits(num_bits_3args(error,
                              static_cast<double>(number_of_elements),
                              static_cast<double>(_num_hash_functions)))
    , data_vec(std::vector<BloomFilterU64<thread_safe>>((num_bits + 63) / 64))
  { }

  [[nodiscard]] size_t num_bits_get(void) const
  {
    return num_bits;
  }

  // returns true if it is already inserted
  bool insert(uint64_t value)
  {
    bool contained = true;
    for (size_t idx = 0; idx < num_hash_functions; idx++)
    {
      const uint64_t hash_value = hash_function(value, idx);
      contained &= set_bit(hash_value % num_bits);
    }
    return contained;
  }

  [[nodiscard]] bool contains(uint64_t value) const
  {
    for (size_t idx = 0; idx < num_hash_functions; idx++)
    {
      const uint64_t hash_value = hash_function(value, idx);
      if (not get_bit(hash_value % num_bits))
      {
        return false;
      }
    }
    return true;
  }

  [[nodiscard]] size_t size_in_bytes(void) const
  {
    return data_vec.size() * sizeof(BloomFilterU64<thread_safe>);
  }

  [[nodiscard]] size_t num_hash_functions_get(void) const
  {
    return num_hash_functions;
  }

  void stats(void) const
  {
    printf("Bloom filter size: %zu bytes\n", size_in_bytes());
  }
};

#endif // BLOOM_FILTER_HPP
