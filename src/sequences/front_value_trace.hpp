#ifndef FRONTVALUE_TRACE_HPP
#define FRONTVALUE_TRACE_HPP
#include <cstddef>
#include <cstdint>
#include <cstdbool>
#include <cassert>
#include <string>

class Backreference
{
  static constexpr const uint8_t ft_eop_mismatch = uint8_t(1),
                                 ft_eop_insertion = (uint8_t(1) << 1),
                                 ft_eop_deletion = (uint8_t(1) << 2);

  uint8_t backbits;
  uint32_t local_matchcount;
  public:
  Backreference(void)
    : backbits(0)
    {}

  Backreference(size_t _local_matchcount)
    : backbits(0)
    , local_matchcount(_local_matchcount)
    {}

  void local_matchcount_set(size_t _local_matchcount)
  {
    assert(_local_matchcount <= UINT32_MAX);
    local_matchcount = static_cast<uint32_t>(_local_matchcount);
  }
  uint32_t local_matchcount_get(void) const noexcept { return local_matchcount;}

  void deletion_set(void) { backbits = ft_eop_deletion; }
  void deletion_add(void) { backbits |= ft_eop_deletion; }
  void insertion_set(void) { backbits = ft_eop_insertion; }
  void insertion_add(void) { backbits |= ft_eop_insertion; }
  void mismatch_set(void) { backbits = ft_eop_mismatch; }
  void mismatch_add(void) { backbits |= ft_eop_mismatch; }

  bool has_deletion(void) const noexcept { return backbits & ft_eop_deletion;}
  bool has_insertion(void) const noexcept { return backbits & ft_eop_insertion;}
  bool has_mismatch(void) const noexcept { return backbits & ft_eop_mismatch;}
  std::string to_string(void) const noexcept
  {
    std::string s{};
    if (has_deletion())
    {
      s += "D";
    }
    if (has_insertion())
    {
      s += "I";
    }
    if (has_mismatch())
    {
      s += "X";
    }
    s += std::to_string(local_matchcount);
    return s;
  }
};

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
  std::string to_string(void) const noexcept
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

  Backreference backreference_get(void) const noexcept
  {
    return backreference;
  }
};
#endif
