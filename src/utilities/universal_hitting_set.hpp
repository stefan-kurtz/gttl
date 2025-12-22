#ifndef UHS_HPP
#define UHS_HPP

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <unordered_set>
#include "utilities/one_hashing_blocked_bloom_filter.hpp"

struct UHSMinimizerKey
{
  bool in_uhs;
  uint64_t hash;

  [[nodiscard]] bool operator < (const UHSMinimizerKey& other)
    const noexcept
  {
    if (in_uhs != other.in_uhs) return in_uhs;
    return hash < other.hash;
  }
};

template <class Hash, class UHS>
struct UHSKeyBuilder
{
  const UHS& uhs;
  static constexpr Hash hash{};

  [[nodiscard]] UHSMinimizerKey operator () (std::string_view kmer)
    const noexcept
  {
    const uint64_t h = hash(kmer);
    return UHSMinimizerKey(uhs.contains(h), h);
  }
};

template <bool thread_safe>
class GttlUHS
{
  private:
  OneHashingBlockedBloomFilter<thread_safe> bloom;
  std::unordered_set<uint64_t> exact;

  public:
  explicit GttlUHS(size_t bloom_blocks,
                   size_t exact_reserve = 0)
    : bloom(bloom_blocks)
  {
    if (exact_reserve != 0) exact.reserve(exact_reserve);
  }

  inline void insert(uint64_t value)
  {
    bloom.insert(value);
    exact.insert(value);
  }

  [[nodiscard]] inline bool may_contain(uint64_t value) const noexcept
  {
    return bloom.contains(value);
  }

  [[nodiscard]] inline bool contains(uint64_t value) const noexcept
  {
    if (not bloom.contains(value)) return false;
    return exact.contains(value);
  }

  [[nodiscard]] size_t bloom_size() const noexcept
  {
    return bloom.size_in_bytes();
  }

  [[nodiscard]] size_t exact_size() const noexcept
  {
    return exact.size();
  }
};

#endif // UHS_HPP
