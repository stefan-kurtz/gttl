#include <iostream>
#include <cstring>
#include "utilities/gttl_line_iterator.hpp"

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
      inputfiles.push_back(std::string(argv[idx]));
    }
    bool all_empty_files = true;
    try
    {
      GttlLineIterator<buf_size> gttl_li(&inputfiles);
      while (gttl_li.next(&buffer))
      {
        all_empty_files = false;
        std::cout << buffer;
        buffer.clear();
      }
    }
    catch (std::string &msg)
    {
      std::cerr << argv[0] << ": file \"" << argv[2] << "\""
                << msg << std::endl;
      haserr = true;
    }
    if (all_empty_files)
    {
      std::cerr << argv[0] << ": file \"" << argv[2]
                << "\", line 1: corrupted sequence" << std::endl;
      haserr = false;
    }
  } else
  {
    for (int idx = 1; idx < argc; idx++)
    {
      bool this_file_is_empty = true;
      try
      {
        GttlLineIterator<buf_size> gttl_li(argv[idx]);
        while (gttl_li.next(&buffer))
        {
          std::cout << buffer;
          this_file_is_empty = false;
          buffer.clear();
        }
      }
      catch (std::string &msg)
      {
        std::cerr << argv[0] << ": file \"" << argv[1] << "\""
                  << msg << std::endl;
        haserr = true;
        break;
      }
      if (this_file_is_empty)
      {
        std::cerr << argv[0] << ": file \"" << argv[1]
                  << "\", line 1: corrupted sequence" << std::endl;
        haserr = true;
        break;
      }
    }
  }
  return haserr ? EXIT_FAILURE : EXIT_SUCCESS;
}
