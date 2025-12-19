/*
** Developed by Henning Lindemann
 */

#ifndef ONE_HASHING_BLOCKED_BLOOM_FILTER_HPP
#define ONE_HASHING_BLOCKED_BLOOM_FILTER_HPP

#include "utilities/bloom_filter_u64.hpp"
#include "utilities/bloom_filter_hash_function.hpp"
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

constexpr size_t OneHashingBlockedBloomFilterBlockSize = 512;
static constexpr std::array<uint64_t, 3> partitions{151, 179, 181};

template <bool ThreadSafe>
class OneHashingBlockedBloomFilter
{
 private:
  class Block
  {
   private:
    std::array<BloomFilterU64<ThreadSafe>, 8> data;

    bool set_bit(uint64_t index)
    {
      return data[index / 64].set_bit(index % 64);
    }

    [[nodiscard]] bool get_bit(uint64_t index) const
    {
      return data[index / 64].get_bit(index % 64);
    }

   public:
    bool insert(uint64_t h)
    {
      bool contained = true;
      uint64_t offset = 0;
      for (const uint64_t p : partitions)
      {
        contained &= set_bit((h % p) + offset);

        offset += p;
      }
      return contained;
    }

    [[nodiscard]] bool contains(uint64_t h) const
    {
      uint64_t offset = 0;
      for (const uint64_t p : partitions)
      {
        if (not get_bit((h % p) + offset))
        {
          return false;
        }

        offset += p;
      }
      return true;
    }
  };
  std::vector<Block> blocks;

 public:
  explicit OneHashingBlockedBloomFilter(uint64_t num_blocks)
    : blocks(std::vector<Block>(num_blocks))
  { }

  bool insert(uint64_t value)
  {
    const uint64_t hash_value = hash_function(value, 0);
    return blocks[hash_value % blocks.size()].insert(hash_value);
  }

  [[nodiscard]] bool contains(uint64_t value) const
  {
    const uint64_t hash_value = hash_function(value, 0);
    return blocks[hash_value % blocks.size()].contains(hash_value);
  }
  [[nodiscard]] size_t size_in_bytes(void) const
  {
    return blocks.size() * 64;
  }
  [[nodiscard]] size_t num_hash_functions_get(void) const
  {
    return 1;
  }
};

#endif // ONE_HASHING_BLOCKED_BLOOM_FILTER_HPP
