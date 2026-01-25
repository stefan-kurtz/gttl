#ifndef BITVECTOR_HPP
#define BITVECTOR_HPP
#include <atomic>
#include <bit>
#include <cassert>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <string>
#include <type_traits>

template <bool thread_safe = false>
class Bitvector
{
  static constexpr const size_t bits = sizeof(uint64_t) * CHAR_BIT;
  static constexpr const uint64_t first_bit = uint64_t(1) << (bits - 1);
  std::conditional_t<thread_safe, std::atomic<uint64_t>, uint64_t> value;

  [[nodiscard]] uint64_t the_bit(size_t idx) const noexcept
  {
    assert(idx < bits);
    return first_bit >> idx;
  }
  public:
  explicit Bitvector(void)
    : value(0)
    {}
  explicit Bitvector(uint64_t _value)
    : value(_value)
    {}
  bool set(size_t idx) requires (not thread_safe)
  {
    assert(idx < bits);
    const bool old = static_cast<bool>(value & the_bit(idx));
    value |= the_bit(idx);
    return old;
  }
  bool set(size_t idx) requires (thread_safe)
  {
    assert(idx < bits);
    return static_cast<bool>(
      value.fetch_or(the_bit(idx), std::memory_order_seq_cst) & the_bit(idx)
    );
  }
  void reset(size_t idx) requires (not thread_safe)
  {
    assert(idx < bits);
    value &= ~the_bit(idx);
  }
  void reset(size_t idx) requires (thread_safe)
  {
    assert(idx < bits);
    value.fetch_and(~the_bit(idx), std::memory_order_seq_cst);
  }
  [[nodiscard]]
  bool operator[](size_t idx) const noexcept requires (not thread_safe)
  {
    assert(idx < bits);
    return static_cast<bool>(value & the_bit(idx));
  }
  [[nodiscard]]
  bool operator[](size_t idx) const noexcept requires (thread_safe)
  {
    assert(idx < bits);
    return static_cast<bool>(
      value.load(std::memory_order_relaxed) & the_bit(idx));
  }
  void operator |= (const Bitvector<false> &other)
    noexcept requires (not thread_safe)
  {
    value |= other.value;
  }
  void operator |= (const Bitvector<true> &other)
    noexcept requires (thread_safe)
  {
    value.fetch_or(other.value.load(std::memory_order_relaxed),
                   std::memory_order_seq_cst);
  }
  bool operator == (const Bitvector<false> &rhs)
    const noexcept requires (not thread_safe)
  {
    return value == rhs.value;
  }
  bool operator == (const Bitvector<true> &rhs)
    const noexcept requires (thread_safe)
  {
    return value.load(std::memory_order_relaxed)
           ==
           rhs.value.load(std::memory_order_relaxed);
  }
  bool operator != (const Bitvector &rhs)
  {
    return !(*this == rhs);
  }
  [[nodiscard]] size_t count(void) const noexcept requires (not thread_safe)
  {
    return std::popcount(value);
  }
  [[nodiscard]] size_t count(void) const noexcept requires (thread_safe)
  {
    return std::popcount(value.load(std::memory_order_relaxed));
  }
  [[nodiscard]] std::string to_string(void)
    const noexcept requires (not thread_safe)
  {
    std::string s{};
    s.reserve(bits);
    for (size_t idx = 0; idx < bits; idx++)
    {
      s += (*this)[idx] ? '1' : '0';
    }
    return s;
  }
  [[nodiscard]] std::string to_string(void)
    const noexcept requires (thread_safe)
  {
    const uint64_t snapshot = value.load(std::memory_order_relaxed);
    std::string s{};
    s.reserve(bits);
    for (size_t idx = 0; idx < bits; ++idx)
    {
      s += (snapshot & the_bit(idx)) ? '1' : '0';
    }
    return s;
  }

  bool set_bit(size_t idx) { return set(idx); }
  bool get_bit(size_t idx) const { return (*this)[idx]; }
};
#endif
