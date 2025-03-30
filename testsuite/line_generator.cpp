#include "utilities/gttl_line_generator.hpp"
#include <cstring>
#include <iostream>

int main(int argc,char *argv[])
{
  // We use 2^20 here, because one of our tests creates a line of size ~985k
  constexpr const int buf_size = 1U << 20U;
  bool haserr = false;

  for (int idx = 1; idx < argc; idx++)
  {
    bool this_file_is_empty = true;
    try
    {
      GttlLineGenerator<buf_size, false> gttl_li(argv[idx]);
      for (const auto& line : gttl_li)
      {
        std::cout << line << '\n';
        this_file_is_empty = false;
      }
    }
    catch (std::string &msg)
    {
      std::cerr << argv[0] << ": file \"" << argv[1] << "\""
                << msg << '\n';
      haserr = true;
      break;
    }
    if (this_file_is_empty)
    {
      std::cerr << argv[0] << ": file \"" << argv[1]
                << "\", line 1: corrupted sequence" << '\n';
      haserr = true;
      break;
    }
  }
  return haserr ? EXIT_FAILURE : EXIT_SUCCESS;
}
