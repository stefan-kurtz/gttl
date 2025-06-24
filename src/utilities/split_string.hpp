#ifndef SPLIT_STRING_HPP
#define SPLIT_STRING_HPP
#include <cstddef>
#include <functional>
#include <string>
#include <cassert>
#include <algorithm>
#include <vector>

template<typename T,T convert(size_t, const std::string &)>
static inline std::vector<T> gttl_split_string(const std::string &str,char sep,
                                               int skip = 1)
{
  assert(skip >= 1);
  auto previous = str.cbegin();
  std::vector<T> result{};
  while (true)
  {
    auto next = std::find(previous, str.cend(),sep);
    if (next < str.cend())
    {
      assert(*next == sep);
      std::string this_string = std::string(previous,next);
      result.push_back(convert(result.size(),this_string));
      if (sep == ' ')
      {
        while (next + 1 < str.cend() and *(next+1) == ' ')
        {
          ++next;
        }
      }
    } else
    {
      std::string this_string = std::string(previous, *(next-1) == '\n'
                                                        ? (next - 1)
                                                        : next);
      result.push_back(convert(result.size(),this_string));
      break;
    }
    previous = next + skip;
  }
  return result;
}

/* This is the same function as the previous except for the fact
   the strings are put into a vector without applying a conversion
   function */

static inline std::string split_string_identity([[maybe_unused]]
                                                 size_t column_idx,
                                                 const std::string &s)
{
  return std::ref(s);
}

static inline std::vector<std::string> gttl_split_string(const std::string &str,
                                                         char sep,
                                                         int skip = 1)
{
  return gttl_split_string<std::string,split_string_identity>(str,sep,skip);
}
#endif
