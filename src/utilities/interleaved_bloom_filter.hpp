/*
** Developed by Henning Lindemann
 */

#ifndef INTERLEAVED_BLOOM_FILTER_HPP
#define INTERLEAVED_BLOOM_FILTER_HPP

#include "sequences/multiseq_factory.hpp"
#include "utilities/bloom_filter.hpp"
#include "utilities/bloom_filter_hash_function.hpp"
#include "utilities/multibitvector.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <numbers>
#include <vector>

template <class T>
class SingleInterleavedBloomFilter
{
  // This gives it access to the private set_bit function
  friend class InterleavedBloomFilter2;
 private:
  const size_t num_hash_functions;
  std::vector<T> data;

  void set_bit(size_t bin, size_t pos)
  {
    data[pos] |= (static_cast<T>(1)) << bin;
  }
 public:
  SingleInterleavedBloomFilter(uint64_t num_bits, size_t _num_hash_functions)
    : num_hash_functions(_num_hash_functions)
    , data(std::vector<T>(num_bits))
  { }


  void insert(size_t bin, uint64_t value)
  {
    for (size_t idx = 0; idx < num_hash_functions; idx++)
    {
      const uint64_t hash_value = hash_function(value, idx);
      // printf("before: %064lb %zu\n",
      //        data[hash_value % data.size()], hash_value %
      //        data.size());
      data[hash_value % data.size()] |= (static_cast<T>(1) << bin);
      // printf("inserted: %064lb %zu\n", data[hash_value % data.size()],
      //         index % data.size());
    }
  }

  [[nodiscard]] T contains(uint64_t value) const
  {
    const uint64_t first_hash_value = hash_function(value, 0);
    T result = data[first_hash_value % data.size()];

    // printf("result: %064lb\n", result);
    for (size_t idx = 1; idx < num_hash_functions; idx++)
    {
      const uint64_t hash_value = hash_function(value, idx);
      result &= data[hash_value % data.size()];
      // printf("result: %064lb\n", result);
    }
    return result;
  }

  [[nodiscard]] size_t size(void) const
  {
    return data.size() * sizeof(T);
  }

  [[nodiscard]] size_t num_bits_get() const
  {
    return data.size();
  }

  [[nodiscard]] size_t num_hash_functions_get(void) const
  {
    return num_hash_functions;
  }
};

class InterleavedBloomFilter
{
 private:
  std::vector<SingleInterleavedBloomFilter<uint64_t>> bloom_filter;
  size_t bins;

 public:
  InterleavedBloomFilter(size_t bins, size_t num_bits,
                         size_t num_hash_functions)
      : bins(bins)
  {
    for (size_t idx = 0; idx < (bins + 63) / 64; idx++)
    {
      bloom_filter.emplace_back(num_bits, num_hash_functions);
    }
  }
  InterleavedBloomFilter(const GttlMultiseqFactory* multiseq_factory,
                         double error)
      : bins(multiseq_factory->size())
  {
    for (size_t idx = 0; idx < (bins + 63) / 64; idx++)
    {
      size_t num_bits = 0;
      size_t num_hash_functions = 0;
      for (size_t j = 0; j < std::min(size_t{64}, bins - idx * 64); j++)
      {
        const size_t number_of_elements = multiseq_factory->at(idx * 64 + j)
                                            ->sequences_total_length_get();
        const double ln_error = std::log(error);

        num_hash_functions =
          std::max(num_hash_functions,
                   static_cast<size_t>(ln_error * -std::numbers::log2e));
        num_bits = std::max<size_t>(num_bits,
                                    static_cast<double>(number_of_elements)
                                      * ln_error * -2.081368);
      }
      bloom_filter.emplace_back(num_bits, num_hash_functions);
    }
  }

  void insert(size_t bin, uint64_t value)
  {
    bloom_filter.at(bin / 64).insert(bin % 64, value);
  }

  Multibitvector<true> contains(uint64_t value)
  {
    Multibitvector<true> result(bins);
    for (size_t idx = 0; idx < (bins + 63) / 64; idx++)
    {
      const uint64_t local_result = bloom_filter.at(idx).contains(value);
      // printf("local_result: %064lb\n", local_result);
      for (size_t j = 0; j < std::min<size_t>(64, idx * 64 - bins); j++)
      {
        if (((local_result >> j) & 1) == 1)
        {
          result.set(idx * 64 + j);
        }
      }
    }
    return result;
  }

  [[nodiscard]] size_t size(void) const
  {
    size_t sum = 0;
    for (const auto & idx : bloom_filter)
    {
      sum += idx.size();
    }
    return sum;
  }

  [[nodiscard]] size_t individual_size(void) const
  {
    size_t max = 0;
    for (const auto & idx : bloom_filter)
    {
      max = std::max(max, idx.size());
    }
    return max;
  }
};

class InterleavedBloomFilter2
{
 private:
  std::vector<SingleInterleavedBloomFilter<uint64_t>> bloom_filter_vec;
  size_t bins;
  BloomFilter<false> insert_bloom_filter;
  bool insert_bloom_filter_is_valid = false;
  size_t current_insert = 0;

 public:
  InterleavedBloomFilter2(GttlMultiseqFactory* multiseq_factory, double error)
      : bins(multiseq_factory->size()),
        insert_bloom_filter(BloomFilter<false>(size_t{0}, size_t{0}))
  {
    for (size_t idx = 0; idx < (bins + 63) / 64; idx++)
    {
      size_t num_bits = 0;
      size_t num_hash_functions = 0;
      for (size_t j = 0; j < std::min(size_t{64}, bins - idx * 64); j++)
      {
        const size_t number_of_elements = multiseq_factory->at(idx * 64 + j)
                                            ->sequences_total_length_get();
        const double ln_error = std::log(error);

        num_hash_functions =
            std::max<size_t>(num_hash_functions,
                             ln_error * -std::numbers::log2e);
        num_bits = std::max<size_t>(
            num_bits,
            (static_cast<double>(number_of_elements)) * ln_error * -2.081368);
      }
      bloom_filter_vec.emplace_back(num_bits, num_hash_functions);
    }
  }

  void insert(size_t bin, uint64_t value)
  {
    if (insert_bloom_filter_is_valid && bin != current_insert)
    {
      flush_insert_bloom_filter();
    }
    if (!insert_bloom_filter_is_valid)
    {
      const size_t num_bits = bloom_filter_vec.at(bin / 64).num_bits_get();
      const size_t num_hash_functions
        = bloom_filter_vec.at(bin / 64).num_hash_functions_get();
      insert_bloom_filter = BloomFilter<false>(num_bits, num_hash_functions);
      current_insert = bin;
      insert_bloom_filter_is_valid = true;
    }
    insert_bloom_filter.insert(value);
  }

  Multibitvector<true> contains(uint64_t value)
  {
    if (insert_bloom_filter_is_valid)
    {
      flush_insert_bloom_filter();
    }
    Multibitvector<true> result(bins);
    for (size_t idx = 0; idx < (bins + 63) / 64; idx++)
    {
      const uint64_t local_result = bloom_filter_vec.at(idx).contains(value);
      // printf("local_result: %064lb\n", local_result);
      for (size_t j = 0; j < std::min<size_t>(64, idx * 64 - bins); j++)
      {
        if (((local_result >> j) & 1) == 1)
        {
          result.set(idx * 64 + j);
        }
      }
    }
    return result;
  }

  [[nodiscard]] size_t size() const
  {
    size_t sum = 0;
    for (const auto & idx : bloom_filter_vec)
    {
      sum += idx.size();
    }
    return sum;
  }

  [[nodiscard]] size_t individual_size() const
  {
    size_t max = 0;
    for (const auto & idx : bloom_filter_vec)
    {
      max = std::max(max, idx.size());
    }
    return max;
  }

  void flush_insert_bloom_filter(void)
  {
    for (size_t idx = 0; idx < insert_bloom_filter.num_bits_get(); idx++)
    {
      if (insert_bloom_filter.get_bit(idx))
      {
        bloom_filter_vec.at(current_insert / 64)
                        .set_bit(current_insert % 64, idx);
      }
    }
    insert_bloom_filter_is_valid = false;
  }
};

#endif // INTERLEAVED_BLOOM_FILTER_HPP
