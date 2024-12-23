#ifndef IS_IN_PATH_HPP
#define IS_IN_PATH_HPP
#include "utilities/split_string.hpp"
#include "utilities/file_exists.hpp"

static inline bool gttl_is_in_PATH(const std::string &prog)
{
  const char *env_p = std::getenv("PATH");
  const auto path_list = gttl_split_string(std::string(env_p),':',1);
  for (auto &path : path_list)
  {
    if (gttl_file_exists(path + "/" + prog))
    {
      return true;
    }
  }
  return false;
}
#endif
