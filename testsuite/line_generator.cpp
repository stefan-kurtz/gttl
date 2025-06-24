#include "utilities/gttl_line_generator.hpp"
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

int main(int argc,char *argv[])
{
  constexpr const int buf_size = 1 << 14;
  std::string buffer{};
  bool haserr = false;

  if (argc > 2 && strcmp(argv[1],"--all") == 0)
  {
    std::vector<std::string> inputfiles{};
    for (int idx = 2; idx < argc; idx++)
    {
      inputfiles.emplace_back(argv[idx]);
    }
    bool all_empty_files = true;
    try
    {
      GttlLineGenerator<buf_size> gttl_li(&inputfiles);
      for (const auto& line : gttl_li)
      {
        all_empty_files = false;
        std::cout << line << '\n';
        buffer.clear();
      }
    }
    catch (const std::exception &msg)
    {
      std::cerr << argv[0] << ": file \"" << argv[2] << "\""
                << msg.what() << '\n';
      haserr = true;
    }
    if (!haserr && all_empty_files)
    {
      std::cerr << argv[0] << ": file \"" << argv[2]
                << "\", line 1: corrupted sequence\n";
      haserr = false;
    }
  } else
  {
    for (int idx = 1; idx < argc; idx++)
    {
      bool this_file_is_empty = true;
      try
      {
        GttlLineGenerator<buf_size> gttl_li(argv[idx]);
        for (const auto& line : gttl_li)
        {
          std::cout << line << '\n';
          this_file_is_empty = false;
          buffer.clear();
        }
      }
      catch (const std::exception &msg)
      {
        std::cerr << argv[0] << ": file \"" << argv[1] << "\""
                  << msg.what() << '\n';
        haserr = true;
        break;
      }
      if (this_file_is_empty)
      {
        std::cerr << argv[0] << ": file \"" << argv[1]
                  << "\", line 1: corrupted sequence\n";
        haserr = true;
        break;
      }
    }
  }
  return haserr ? EXIT_FAILURE : EXIT_SUCCESS;
}
