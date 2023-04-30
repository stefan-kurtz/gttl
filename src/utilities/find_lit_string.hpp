#ifndef FIND_LIT_STRING_HPP
#define FIND_LIT_STRING_HPP
#include <array>
#include <cstddef>

static inline constexpr bool gttl_compare_string_literals(size_t idx,
                                                          const char *s,
                                                          const char *t)
{
  return (s[idx] == '\0' && t[idx] == '\0')
           ? true
           : (s[idx] == '\0' or t[idx] == '\0') ? false
                                                : s[idx] == t[idx] &&
                                                  gttl_compare_string_literals(
                                                    idx+1,s,t);
}

static inline constexpr bool gttl_compare_string_literals(const char *s,
                                                          const char *t)
{
  return gttl_compare_string_literals(0,s,t);
}

template<size_t num_values>
using GttlLitStringArray = std::array<const char *,num_values>;

template<size_t num_values>
static inline constexpr size_t find_lit_string_at_compile_time(
  size_t idx,
  const GttlLitStringArray<num_values> &arr,
  const char *v)
{
  return idx == num_values
           ? num_values
           : (gttl_compare_string_literals(v,arr[idx])
                ? idx
                : find_lit_string_at_compile_time<num_values>(idx+1,arr,v));
}

template<size_t num_values>
static inline constexpr size_t find_lit_string_at_compile_time(
  const GttlLitStringArray<num_values> &arr,
  const char *v)
{
  return find_lit_string_at_compile_time<num_values>(0,arr,v);
}
#endif
