#ifndef HAS_SUFFIX_OR_PREFIX_HPP
#define HAS_SUFFIX_OR_PREFIX_HPP
#include <string>
#include <vector>
static inline bool gttl_has_prefix(const std::string &s,
                                   const std::string &prefix)
{
  return s.size() >= prefix.size() and s.substr(0,prefix.size()) == prefix;
}

static inline bool gttl_has_suffix(const std::string &s,
                                   const std::string &suffix)
{
  return s.size() >= suffix.size() and
         s.substr(s.size() - suffix.size()) == suffix;
}

static inline bool gttl_has_suffix_with_extension(const std::string &to_check,
                                                  const std::string &suffix,
                                                  const std::string &extension)
{
  return gttl_has_suffix(to_check, suffix) or
         gttl_has_suffix(to_check, suffix + extension);
}

static inline bool gttl_has_any_suffix_with_extension(
  const std::string &to_check,
  const std::vector<std::string> &suffixes,
  const std::string &extension)
{
  for (auto &&suffix : suffixes)
  {
    if (gttl_has_suffix_with_extension(to_check, suffix, extension))
    {
      return true;
    }
  }
  return false;
}
#endif
