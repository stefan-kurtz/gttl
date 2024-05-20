#ifndef HAS_SUFFIX_OR_PREFIX_HPP
#define HAS_SUFFIX_OR_PREFIX_HPP
#include <string>
bool inline gttl_has_prefix(const std::string &s,const std::string &prefix)
{
  return s.size() >= prefix.size() and s.substr(0,prefix.size()) == prefix;
}

bool inline gttl_has_suffix(const std::string &s,const std::string &suffix)
{
  return s.size() >= suffix.size() and 
         s.substr(s.size() - suffix.size()) == suffix;
}
#endif
