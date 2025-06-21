#ifndef FRONT_VALUE_MATCH_HISTORY_HPP
#define FRONT_VALUE_MATCH_HISTORY_HPP
#include <cstddef>
#include <cstdint>
#include <cassert>
#include <string>
#include "utilities/bit_sequence.hpp"
#include "utilities/mathsupport.hpp"
#include "utilities/unused.hpp"

class FrontValueMatchHistory
{
 private:
  size_t row;
  uint64_t mismatch_history;

 public:
  static constexpr const bool unit_cost_model = true;
  FrontValueMatchHistory(void) : row(0), mismatch_history(~uint64_t(0)) {}
  FrontValueMatchHistory(size_t _row)
    : row(_row)
    , mismatch_history(uint64_t(0))
  {}
  uint64_t match_history_get(void) const noexcept { return ~mismatch_history; }
  bool operator>=(size_t scalar) const noexcept { return row >= scalar; }
  bool operator<=(const FrontValueMatchHistory &fh) const noexcept
  {
    return row <= fh.row;
  }
  bool operator<(size_t scalar) const noexcept { return row < scalar; }
  void operator+=(size_t match_length)
  {
#ifdef SKDEBUG
    std::cout << bit_sequence2string<uint64_t>(~mismatch_history)
              << "+=(match_length=" + std::to_string(match_length) << ") gives"
              << std::endl;
#endif
    row += match_length;
    mismatch_history = (match_length >= size_t(64))
                           ? uint64_t(0)
                           : (mismatch_history << match_length);
#ifdef SKDEBUG
    std::cout << bit_sequence2string<uint64_t>(~mismatch_history) << std::endl;
#endif
  }
  void operator+=(int value)
  {
#ifdef SKDEBUG
    std::cout << bit_sequence2string<uint64_t>(~mismatch_history)
              << "+=(value=" + std::to_string(value) << ") gives"
              << std::endl;
#endif
    row += value;
    mismatch_history = (mismatch_history << 1) | uint64_t(1);
#ifdef SKDEBUG
    std::cout << bit_sequence2string<uint64_t>(~mismatch_history) << std::endl;
#endif
  }
  size_t operator+(size_t ell) const noexcept { return row + ell; }
  std::string to_string(void) const noexcept
  {
    std::string cs{};
    cs += std::to_string(row);
    const uint64_t match_history = match_history_get();
    cs += "\t" + bit_sequence2string<uint64_t>(match_history);
    cs += "\t" + std::to_string(popcount_uint64_t(match_history));
    return cs;
  }
  size_t row_get(size_t ulen) const noexcept { return std::min(row, ulen); }
  size_t aligned_len_get(int32_t diag_idx,size_t ulen,
                         GTTL_UNUSED size_t vlen) const noexcept
  {
    const size_t dest_row = row_get(ulen);
    assert(diag_idx >= 0 || static_cast<size_t>(-diag_idx) <= 2 * dest_row);
    return 2 * dest_row + diag_idx;
  }
};
#endif
