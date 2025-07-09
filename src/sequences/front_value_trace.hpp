#ifndef FRONT_VALUE_TRACE_HPP
#define FRONT_VALUE_TRACE_HPP
#include <cstddef>
#include <string>
#include "sequences/backreference.hpp"

class FrontValueTrace
{
  private:
  size_t row;
  Backreference backreference;

  public:
  FrontValueTrace(void)
    : row(0)
    , backreference({})
    {}

  FrontValueTrace(size_t _row)
    : row(_row)
    , backreference(_row)
    {}

  bool operator>=(size_t scalar) const noexcept { return row >= scalar; }
  bool operator==(size_t scalar) const noexcept { return row == scalar; }
  bool operator<=(const FrontValueTrace &fh) const noexcept
  {
    return row <= fh.row;
  }
  bool operator<(size_t scalar) const noexcept { return row < scalar; }
  void operator+=(size_t match_length)
  {
    row += match_length;
    backreference.local_matchcount_set(match_length);
  }
  void operator+=(int value) { row += value; }
  size_t operator+(size_t ell) const noexcept { return row + ell; }
  [[nodiscard]] std::string to_string(void) const noexcept
  {
    std::string cs{};
    cs += std::to_string(row);
    return cs;
  }
  void deletion_set(void) { backreference.deletion_set(); }
  void deletion_add(void) { backreference.deletion_add(); }
  void insertion_set(void) { backreference.insertion_set(); }
  void insertion_add(void) { backreference.insertion_add(); }
  void mismatch_set(void) { backreference.mismatch_set(); }
  void mismatch_add(void) { backreference.mismatch_add(); }

  [[nodiscard]] Backreference backreference_get(void) const noexcept
  {
    return backreference;
  }
};
#endif
