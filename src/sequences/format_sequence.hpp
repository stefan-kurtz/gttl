#ifndef FORMAT_SEQUENCE_HPP
#define FORMAT_SEQUENCE_HPP
#include <iostream>
#include <string_view>
#include <cstddef>

static inline void gttl_format_sequence(const std::string_view &sequence,
                                        size_t line_width)
{
  if (line_width == 0)
  {
    std::cout << sequence << '\n';
  } else
  {
    size_t remain = sequence.size();
    size_t offset = 0;
    while (true)
    {
      if (line_width < remain)
      {
        std::cout << sequence.substr(offset, line_width) << '\n';
        remain -= line_width;
        offset += line_width;
      } else
      {
        std::cout << sequence.substr(offset, remain) << '\n';
        break;
      }
    }
  }
}
#endif
