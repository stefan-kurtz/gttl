#include "utilities/gttl_binary_read.hpp"
#include <iostream>
#include <stdexcept>

int main(int argc, char *argv[])
{
  for (int idx = 1; idx < argc; idx++)
  {
    try
    {
      BinaryFileReader<char> binary_file_reader(argv[1]);
      for (const auto &entry : binary_file_reader)
      {
        std::cout << entry;
      }
    }
    catch (const std::out_of_range &err)
    {
      std::cerr << argv[0] << ": file \"" << argv[idx] << "\": " << err.what()
                << '\n';
      return EXIT_FAILURE;
    }
    catch (const std::runtime_error &err)
    {
      std::cerr << argv[0] << ": file \"" << argv[idx] << "\": " << err.what()
                << '\n';
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}
