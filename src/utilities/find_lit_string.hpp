#ifndef FIND_LIT_STRING_HPP
#define FIND_LIT_STRING_HPP
#include <initializer_list>
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

using GttlLitStringInitializerList = std::initializer_list<const char *>;

static inline consteval size_t gttl_find_lit_string_at_compile_time(
  size_t idx,
  const GttlLitStringInitializerList &arr,
  const char *v)
{
  return idx == arr.size()
           ? idx
           : (gttl_compare_string_literals(v,*(arr.begin() + idx))
                ? idx
                : gttl_find_lit_string_at_compile_time(idx+1,arr,v));
}

static inline consteval size_t gttl_find_lit_string_at_compile_time(
  const GttlLitStringInitializerList &arr,
  const char *v)
{
  return gttl_find_lit_string_at_compile_time(0,arr,v);
}
#endif
