#include <cassert>
#include <cstdlib>
#include <iostream>
#include "utilities/string_of_digits.hpp"

#define DIGITS "1234567890"
#define DIGITS2 DIGITS DIGITS
#define DIGITS4 DIGITS2 DIGITS2
#define DIGITS8 DIGITS4 DIGITS4
#define DIGITS16 DIGITS8 DIGITS8
#define DIGITS32 DIGITS16 DIGITS16
#define DIGITS64 DIGITS32 DIGITS32

constexpr const char *is_digits[] = {"0",      "12345",  DIGITS,
                                     DIGITS2,  DIGITS4,  DIGITS8,
                                     DIGITS16, DIGITS32, DIGITS64};
constexpr const char *is_not_digits[] = {"",       "a",      "12345a", "a12345",
                                         "123a45", "123!45", "-12345", "765.98",
                                         "818 12", " 12345 "};

int main()
{
  for (const auto str : is_digits)
  {
    if (!string_of_digits(str))
    {
      std::cerr << '"' << str << '"'
                << " was not recognized as a string of digits!\n";
      return EXIT_FAILURE;
    }
  }
  for (const auto str : is_not_digits)
  {
    if (string_of_digits(str))
    {
      std::cerr << '"' << str << '"'
                << " was recognized as a string of digits!\n";
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}
