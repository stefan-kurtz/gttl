#ifndef SPLIT_STRING_HPP
#define SPLIT_STRING_HPP
#include <string>
#include <cstdbool>

template<typename T,T convert(size_t,const std::string &)>
static inline std::vector<T> split_string(const std::string &str,char sep,
                                          int skip = 1)
{
  assert(skip >= 1);
  auto previous = str.begin();
  std::vector<T> result{};
  while (true)
  {
    auto next = std::find(previous, str.end(),sep);
    std::string this_string = std::string(previous,next);
    result.push_back(convert(result.size(),this_string));
    if (next == str.end())
    {
      break;
    }
    previous = next + skip;
  }
  return result;
}

/* This is the same function as the previous except for the fact
   the strings are put into a vector without applying a conversion
   function */

static inline std::vector<std::string> split_string(const std::string &str,
                                                    char sep,
                                                    int skip = 1)
{
  assert(skip >= 1);
  auto previous = str.begin();
  std::vector<std::string> result{};
  while (true)
  {
    auto next = std::find(previous, str.end(),sep);
    std::string this_string = std::string(previous,next);
    result.push_back(this_string);
    if (next == str.end())
    {
      break;
    }
    previous = next + skip;
  }
  return result;
}
#endif
