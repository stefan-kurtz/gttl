#ifndef IS_IN_PATH_HPP
#define IS_IN_PATH_HPP
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <string>
#include "utilities/split_string.hpp"

static inline bool gttl_is_in_PATH(const std::string &prog)
{
  const char *const env_p = std::getenv("PATH");
  const auto path_list = gttl_split_string(std::string(env_p),':',1);
  return std::ranges::any_of(path_list,
                             [&](const auto& path)
                             { return std::filesystem::exists(
                                 std::filesystem::path(path) / prog
                                );
                             });
}
#endif
