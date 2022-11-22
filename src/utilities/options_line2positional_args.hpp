#ifndef OPTIONS_LINE2POSITIONAL_ARGS_HPP
#define OPTIONS_LINE2POSITIONAL_ARGS_HPP

#include <vector>
#include <string>
#include <cassert>
#include "utilities/split_string.hpp"

static inline std::vector<std::string>
  gttl_options_line2positional_args(const std::string &options_line,
                                    const std::string &sep)
{
  constexpr const int skip = 1;
  std::vector<std::string> argv = split_string(options_line,' ',skip);
  size_t idx = 0;
  for (auto &&arg : argv)
  {
    if (arg == sep)
    {
      break;
    }
    idx++;
  }
  assert (idx < argv.size());
  idx++;
  std::vector<std::string> positional_args{};
  while (idx < argv.size())
  {
    positional_args.push_back(argv[idx]);
    idx++;
  }
  return positional_args;
}
#endif
