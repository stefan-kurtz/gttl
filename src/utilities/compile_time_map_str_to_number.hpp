#ifndef COMPILE_TIME_MAP_STR_TO_NUMBER_HPP
#define COMPILE_TIME_MAP_STR_TO_NUMBER_HPP

#include "utilities/gttl_string_literal.hpp"
#include "utilities/find_lit_string.hpp"

template<const GttlLitStringInitializerList &string_list>
struct CompileTimeMapStrToNumber
{
  template<GttlStringLiteral key>
  [[nodiscard]] constexpr size_t get(void) const noexcept
  {
    return gttl_find_lit_string_at_compile_time(string_list, key.value);
  }
};
#endif
