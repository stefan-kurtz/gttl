#include "utilities/gttl_binary_read.hpp"
#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>

int main(int argc, char *argv[])
{
  for (int idx = 1; idx < argc; idx++)
  {
    try
    {
      const std::string inputfile{argv[idx]};
      const BinaryFileReader<char> binary_file_reader(inputfile);
      for (const auto &entry : binary_file_reader)
      {
        std::cout << entry;
      }
    }
    catch (const std::exception &err)
    {
      std::cerr << argv[0] << ": file \"" << argv[idx] << "\": " << err.what()
                << '\n';
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}
