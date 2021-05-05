#include <filesystem>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <vector>
#include "utilities/str_format.hpp"

template<typename T>
std::vector<T> gttl_read_vector(const char *filename)
{
  // Open the stream to 'lock' the file.
  std::ifstream instream(filename, std::ios::in | std::ios::binary);

  if (instream.fail())
  {
    StrFormat msg("cannot open file %s",filename);
    throw std::ios_base::failure(msg.str());
  }
  const auto size_of_file = std::filesystem::file_size(filename);
  if (size_of_file % sizeof(T) != 0)
  {
    StrFormat msg("file %s contains %lu bytes which is not a multiple of %lu",
                  filename,size_of_file,sizeof(T));
    throw std::ios_base::failure(msg.str());
  }
  size_t num_values = size_of_file/sizeof(T);
  std::vector<T> vec(num_values);
  if (!instream.read(reinterpret_cast<char*>(vec.data()), size_of_file))
  {
    StrFormat msg("cannot only read %lu bytes from file %s",
                  instream.gcount(),filename);
    throw std::ios_base::failure(msg.str());
  }
  return vec;
}
