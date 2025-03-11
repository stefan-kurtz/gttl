#ifndef STRING_OF_DIGITS_HPP
#define STRING_OF_DIGITS_HPP
#include <algorithm>
#include <string>

inline bool string_of_digits(const std::string &s)
{
  return not s.empty() and std::all_of(s.begin(), s.end(), ::isdigit);
}
#endif
