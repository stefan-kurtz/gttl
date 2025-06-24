#ifndef IS_IN_PATH_HPP
#define IS_IN_PATH_HPP
#include <cstdlib>
#include <filesystem>
#include <string>
#include "utilities/split_string.hpp"

static inline bool gttl_is_in_PATH(const std::string &prog)
{
  const char *env_p = std::getenv("PATH");
  const auto path_list = gttl_split_string(std::string(env_p),':',1);
  for (auto &path : path_list)
  {
    if (std::filesystem::exists(path + "/" + prog))
    {
      return true;
    }
  }
  return false;
}
#endif
