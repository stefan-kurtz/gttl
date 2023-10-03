#ifndef SPLIT_STRING_HPP
#define SPLIT_STRING_HPP
#include <string>
#include <cstdbool>
#include <cassert>
#include <algorithm>
#include <vector>

template<typename T,T convert(size_t,const std::string &)>
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
      result.push_back(convert(this_string));
    } else
    {
      assert (next == str.cend() && *(next-1) == '\n');
      std::string this_string = std::string(previous,next - 1);
      result.push_back(convert(this_string));
      break;
    }
    previous = next + skip;
  }
  return result;
}

/* This is the same function as the previous except for the fact
   the strings are put into a vector without applying a conversion
   function */

static inline std::vector<std::string> gttl_split_string(const std::string &str,
                                                         char sep,
                                                         int skip = 1)
{
  assert(skip >= 1);
  auto previous = str.cbegin();
  std::vector<std::string> result{};
  while (true)
  {
    auto next = std::find(previous, str.cend(),sep);
    if (next < str.cend())
    {
      assert(*next == sep);
      std::string this_string = std::string(previous,next);
      result.push_back(this_string);
    } else
    {
      assert(next == str.cend() && *(next-1) == '\n');
      std::string this_string = std::string(previous,next - 1);
      result.push_back(this_string);
      break;
    }
    previous = next + skip;
  }
  return result;
}
#endif
