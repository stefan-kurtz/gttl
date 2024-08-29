#ifndef CONCATENATE_STRINGS_HPP
#define CONCATENATE_STRINGS_HPP

#include <vector>
#include <string>

template<class Iter>
static inline std::string gttl_concatenate_strings(Iter first, Iter last,
                                                   const std::string &sep)
{
  std::string string_from_vec;
  for (auto it = first; it != last; ++it)
  {
    if (string_from_vec.size() > 0)
    {
      string_from_vec += sep;
    }
    string_from_vec += *it;
  }
  return string_from_vec;
}
#endif
