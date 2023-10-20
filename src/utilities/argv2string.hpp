#ifndef ARGV2STRING_HPP
#define ARGV2STRING_HPP
#include <string>

static inline std::string gttl_argv2string(const char *prefix,
                                           int argc, char *argv[])
{
  std::string s(prefix);
  for (int idx = 0; idx < argc; idx++)
  {
    s += " " + std::string(argv[idx]);
  }
  return s;
}
#endif
