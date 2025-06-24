#include <cstdlib>
#include <exception>
#include <string>
#include <iostream>
#include "utilities/gttl_file_open.hpp"

int main(int argc, char* argv[])
{
  if(argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " <filename>\n";
    exit(EXIT_FAILURE);
  }try
  {
    std::string x = gttl_read_file(argv[1]);
    std::cout << x;
  } catch(const std::exception &err)
  {
    std::cerr << argv[0] << ": " << err.what() << '\n';
    exit(EXIT_FAILURE);
  }
}
