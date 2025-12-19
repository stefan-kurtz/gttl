/*
** Developed by Henning Lindemann
 */

#ifndef BLOCKED_BLOOM_FILTER_HPP
#define BLOCKED_BLOOM_FILTER_HPP

#include "utilities/bloom_filter_u64.hpp"
#include "utilities/bloom_filter_hash_function.hpp"
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

constexpr size_t BlockedBloomFilterBlockSize = 512;

template <bool thread_safe>
class BlockedBloomFilter
{
 private:
  class Block
  {
   private:
    std::array<BloomFilterU64<thread_safe>, 8> data;
    bool set_bit(uint64_t index)
    {
      return data[index / 64].set_bit(index % 64);
    }
    [[nodiscard]] bool get_bit(uint64_t index) const
    {
      return data[index / 64].get_bit(index % 64);
    }
   public:
    bool insert(uint64_t value, size_t num_hash_functions)
    {
      bool contained = true;
      for (size_t idx = 0; idx < num_hash_functions; idx++)
      {
        const uint64_t hash_value = hash_function(value, idx);
        contained &= set_bit(hash_value % BlockedBloomFilterBlockSize);
      }
      return contained;
    }

    [[nodiscard]] bool contains(uint64_t value, size_t num_hash_functions) const
    {
      for (size_t idx = 0; idx < num_hash_functions; idx++)
      {
        const uint64_t hash_value = hash_function(value, idx);
        if (!get_bit(hash_value % BlockedBloomFilterBlockSize))
        {
          return false;
        }
      }
      return true;
    }
  };
  std::vector<Block> blocks;
  size_t num_hash_functions;

 public:
  BlockedBloomFilter(uint64_t num_blocks, size_t _num_hash_functions)
    : blocks(std::vector<Block>(num_blocks))
    , num_hash_functions(_num_hash_functions)
  { }

  void insert(uint64_t value)
  {
    const uint64_t hash_value = hash_function(value, 0);

    blocks[hash_value % blocks.size()].insert(value, num_hash_functions);
  }

  [[nodiscard]] bool contains(uint64_t value) const
  {
    const uint64_t hash_value = hash_function(value, 0);
    return blocks[hash_value % blocks.size()].contains(value,
                                                       num_hash_functions);
  }

  [[nodiscard]] size_t size_in_bytes(void) const
  {
    return blocks.size() * 64;
  }

  [[nodiscard]] size_t num_hash_functions_get(void) const
  {
    return num_hash_functions;
  }
};

#endif // BLOCKED_BLOOM_FILTER_HPP
