#ifndef BACKREFERENCE_HPP
#define BACKREFERENCE_HPP
#include <cstddef>
#include <cstdint>
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
  {
    assert(_local_matchcount <= UINT32_MAX);
  }

  void local_matchcount_set(size_t _local_matchcount)
  {
    assert(_local_matchcount <= UINT32_MAX);
    local_matchcount = static_cast<uint32_t>(_local_matchcount);
  }
  [[nodiscard]] uint32_t local_matchcount_get(void) const noexcept
  {
    return local_matchcount;
  }

  void deletion_set(void) { backbits = ft_eop_deletion; }
  void deletion_add(void) { backbits |= ft_eop_deletion; }
  void insertion_set(void) { backbits = ft_eop_insertion; }
  void insertion_add(void) { backbits |= ft_eop_insertion; }
  void mismatch_set(void) { backbits = ft_eop_mismatch; }
  void mismatch_add(void) { backbits |= ft_eop_mismatch; }

  [[nodiscard]] bool has_deletion(void) const noexcept
  {
    return backbits & ft_eop_deletion;
  }
  [[nodiscard]] bool has_insertion(void) const noexcept
  {
    return backbits & ft_eop_insertion;
  }
  [[nodiscard]] bool has_mismatch(void) const noexcept
  {
    return backbits & ft_eop_mismatch;
  }
  [[nodiscard]] std::string to_string(void) const noexcept
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
#endif
