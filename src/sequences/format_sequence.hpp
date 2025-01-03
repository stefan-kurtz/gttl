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
    std::cout << sequence << std::endl;
  } else
  {
    size_t remain = sequence.size(),
           offset = 0;
    while (true)
    {
      if (line_width < remain)
      {
        std::cout << sequence.substr(offset,line_width) << std::endl;
        remain -= line_width;
        offset += line_width;
      } else
      {
        std::cout << sequence.substr(offset,remain) << std::endl;
        break;
      }
    }
  }
}
#endif
