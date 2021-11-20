#include <iostream>
#include <cstring>
#include "utilities/gttl_line_iterator.hpp"

int main(int argc,char *argv[])
{
  constexpr const int buf_size = 1 << 14;

  if (argc > 2 && strcmp(argv[1],"--all") == 0)
  {
    std::vector<std::string> inputfiles{};
    for (int idx = 2; idx < argc; idx++)
    {
      inputfiles.push_back(std::string(argv[idx]));
    }
    GttlLineIterator<buf_size> gttl_li(&inputfiles);
    std::string buffer{};
    while (gttl_li.next(&buffer))
    {
      std::cout << buffer;
      buffer.clear();
    }
  } else
  {
    for (int idx = 1; idx < argc; idx++)
    {
      GttlLineIterator<buf_size> gttl_li(argv[idx]);
      std::string buffer{};
      while (gttl_li.next(&buffer))
      {
        std::cout << buffer;
        buffer.clear();
      }
    }
  }
  return EXIT_SUCCESS;
}
