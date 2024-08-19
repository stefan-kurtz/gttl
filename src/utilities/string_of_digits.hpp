#ifndef STRING_OF_DIGITS_HPP
#define STRING_OF_DIGITS_HPP
#include <string>
#include <cctype>

inline bool string_of_digits(const std::string &s)
{
  return not s.empty() and std::all_of(s.begin(), s.end(),
                                       [](unsigned char cc)
                                       {
                                         return std::isdigit(cc);
                                       }
                                      );
}
#endif
