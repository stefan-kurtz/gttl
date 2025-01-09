#include <cstdlib>
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
  } catch(const std::string &e)
  {
    std::cerr << e << std::endl;
    exit(EXIT_FAILURE);
  }
}
