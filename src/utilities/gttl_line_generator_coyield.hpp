#ifndef GTTL_LINE_GENERATOR_COYIELD_HPP
#define GTTL_LINE_GENERATOR_COYIELD_HPP

#include <fstream>
#include <generator>
#include <string>

std::generator<std::string> gttl_read_lines(const std::string &file_path)
{
  std::ifstream file(file_path);
  if(!file.is_open())
  {
    throw std::ios_base::failure("Failed to open file: " + file_path);
  }
  std::string line;
  while(std::getline(file, line))
  {
    co_yield line;
  }
}

#endif //GTTL_LINE_GENERATOR_COYIELD_HPP
