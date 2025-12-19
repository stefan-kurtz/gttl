/*
** Developed by Henning Lindemann
 */

#ifndef BLOOM_FILTER_U64_HPP
#define BLOOM_FILTER_U64_HPP

#include <atomic>
#include <cassert>
#include <cstdint>

template <bool thread_safe>
class BloomFilterU64
{
};

template <>
class BloomFilterU64<false>
{
  uint64_t data;

  public:
  bool set_bit(uint64_t index)
  {
    assert(index < 64);
    const bool old = static_cast<bool>((data >> index) & uint64_t{1});
    data |= uint64_t{1} << index;
    return old;
  }

  [[nodiscard]] bool get_bit(uint64_t index) const
  {
    assert(index < 64);
    return static_cast<bool>((data >> index) & uint64_t{1});
  }
};

template <>
class BloomFilterU64<true>
{
  std::atomic<uint64_t> data;

  public:
  bool set_bit(uint64_t index)
  {
    assert(index < 64);
    return static_cast<bool>((data.fetch_or(uint64_t{1} << index,
                                            std::memory_order_seq_cst) >> index)
                               & uint64_t{1});
  }

  [[nodiscard]] bool get_bit(uint64_t index) const
  {
    assert(index < 64);
    return static_cast<bool>((data.load(std::memory_order_relaxed) >> index)
                               & uint64_t{1});
  }
};


#endif // BLOOM_FILTER_U64_HPP
