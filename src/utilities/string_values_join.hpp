#ifndef STRING_VALUES_JOIN_HPP
#define STRING_VALUES_JOIN_HPP

#include <string>
#include <type_traits>

template<class FwdIter>
static inline std::string string_values_join(const char *sep,
                                             FwdIter begin,
                                             FwdIter end)
{
  if (begin >= end)
  {
    return std::string("");
  }
  FwdIter it = begin;
  std::string out_string{*it};
  ++it;
  for (; it  != end; ++it)
  {
    auto v = *it;
    using T = decltype(v);
    if constexpr (std::is_same_v<T, std::string> or
                  std::is_same_v<T, const char *>)
    {
      out_string += std::string(sep) + v;
    } else
    {
      out_string += std::string(sep) + std::to_string(v);
    }
  }
  return out_string;
}
#endif
