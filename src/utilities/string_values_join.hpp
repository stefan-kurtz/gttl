#ifndef STRING_VALUES_JOIN_HPP
#define STRING_VALUES_JOIN_HPP

#include <string>
#include <vector>

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
    out_string += std::string(sep) + *it;
  }
  return out_string;
}
#endif
