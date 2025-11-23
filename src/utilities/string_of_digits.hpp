#ifndef STRING_OF_DIGITS_HPP
#define STRING_OF_DIGITS_HPP
#include <algorithm>
#include <cctype>
#include <string>

inline bool string_of_digits(const std::string &s)
{
  return not s.empty() and std::ranges::all_of(s, ::isdigit);
}
#endif
