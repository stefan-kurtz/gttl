#ifndef SPLIT_STRING_HPP
#define SPLIT_STRING_HPP
#include <string>
#include <cstdbool>

template<typename T,T convert(const std::string &)>
static std::vector<T> split_string(const std::string &str, char sep)
{
  auto previous = str.begin();
  std::vector<T> result{};
  while (true)
  {
    auto next = std::find(previous, str.end(),sep);
    if (next == str.end())
    {
      std::string this_string = std::string(previous,next-1);
      result.push_back(convert(this_string));
      break;
    }
    std::string this_string = std::string(previous,next);
    result.push_back(convert(this_string));
    previous = next + 1;
  }
  return result;
}
#endif
