#include "utilities/gttl_binary_read.hpp"
#include <iostream>
#include <string>

int main(int argc, char *argv[])
{
  for (int idx = 1; idx < argc; idx++)
  {
    try
    {
      BinaryFileIterator<char> gttl_bi(argv[1]);
      for (const auto &entry : gttl_bi)
      {
        std::cout << entry;
      }
    }
    catch (std::string &msg)
    {
      std::cerr << argv[0] << ": file \"" << argv[1] << "\"" << msg << '\n';
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}
