#ifndef READ_VECTOR_HPP
#define READ_VECTOR_HPP

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <ios>
#include "utilities/file_size.hpp"
#include "utilities/str_format.hpp"

template<typename T>
std::vector<T> gttl_read_vector(const char *filename)
{
  const size_t size_of_file = gttl_file_size(filename);
  if (size_of_file % sizeof(T) != 0)
  {
    StrFormat msg("file %s contains %lu bytes which is not a multiple of %lu",
                  filename,size_of_file,sizeof(T));
    throw msg.str();
  }
  const size_t num_values = size_of_file/sizeof(T);
  std::vector<T> vec(num_values);
  // Open the stream to 'lock' the file.
  std::ifstream instream(filename, std::ios::in | std::ios::binary);
  if (instream.fail())
  {
    StrFormat msg("cannot open file %s",filename);
    throw msg.str();
  }
  if (!instream.read(reinterpret_cast<char*>(vec.data()), size_of_file))
  {
    StrFormat msg("cannot only read %lu bytes from file %s",
                  instream.gcount(),filename);
    throw msg.str();
  }
  return vec;
}
#endif
