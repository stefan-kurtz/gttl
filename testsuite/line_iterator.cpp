#include <iostream>
#include "utilities/gttl_line_iterator.hpp"

int main(int argc,char *argv[])
{
  for (int idx = 1; idx < argc; idx++)
  {
    constexpr const int buf_size = 1 << 14;
    GttlLineIterator<buf_size> gttl_li(argv[idx]);
    std::string buffer{};
    while (gttl_li.next(&buffer))
    {
      std::cout << buffer;
      buffer.clear();
    }
  }
}
