#ifndef HAS_SUFFIX_OR_PREFIX_HPP
#define HAS_SUFFIX_OR_PREFIX_HPP
#include <string>
#include <vector>
static inline bool gttl_has_suffix_with_extension(const std::string &to_check,
                                                  const std::string &suffix,
                                                  const std::string &extension)
{
  return to_check.ends_with(suffix) or
         to_check.ends_with(suffix + extension);
}

static inline bool gttl_has_any_suffix_with_extension(
  const std::string &to_check,
  const std::vector<std::string> &suffixes,
  const std::string &extension)
{
  for (const std::string& suf : suffixes)
  {
    if (gttl_has_suffix_with_extension(to_check, suf, extension))
    {
      return true;
    }
  }
  return false;
}
#endif
