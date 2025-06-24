#ifndef READ_VECTOR_HPP
#define READ_VECTOR_HPP

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <filesystem>
#include <ios>
#include "utilities/file_size.hpp"
#include "utilities/str_format.hpp"
#include "utilities/gttl_line_generator.hpp"
#include "utilities/has_suffix_or_prefix.hpp"

template<typename T>
std::vector<T> gttl_read_vector(const char *filename)
{
  if (gttl_has_suffix(filename,".gz"))
  {
    std::string file_content;
    GttlLineGenerator line_get(filename);
    for (auto &&line : line_get)
    {
      file_content += line;
      file_content += "\n";
    }
    if (file_content.size() % sizeof(T) != 0)
    {
      const StrFormat msg("file %s contains %zu bytes which is not a "
                          "multiple of %zu",
                          filename,
                          file_content.size(),
                          sizeof(T));
      throw std::ios_base::failure(msg.str());
    }
    std::vector<T> vec(reinterpret_cast<const T*>(file_content.data()),
                       reinterpret_cast<const T*>(file_content.data() +
                                                  file_content.size()));
    return vec;
  }
  if (not std::filesystem::exists(filename))
  {
    throw std::ios_base::failure(std::string("file \"") + filename
                                 + std::string("\" does not exist"));
  }
  const size_t size_of_file = gttl_file_size(filename);
  if (size_of_file % sizeof(T) != 0)
  {
    const StrFormat msg("file %s contains %zu bytes which is not a "
                        "multiple of %zu",
                        filename,
                        size_of_file,
                        sizeof(T));
    throw std::ios_base::failure(msg.str());
  }
  // Open the stream to 'lock' the file.
  std::ifstream instream(filename, std::ios::in | std::ios::binary);
  if (instream.fail())
  {
    const StrFormat msg("cannot open file %s", filename);
    throw std::ios_base::failure(msg.str());
  }
  const size_t num_values = size_of_file/sizeof(T);
  std::vector<T> vec(num_values);
  if (!instream.read(reinterpret_cast<char*>(vec.data()), size_of_file))
  {
    const StrFormat msg("cannot only read %zu bytes from file %s",
                        instream.gcount(),
                        filename);
    throw std::ios_base::failure(msg.str());
  }
  return vec;
}

template<typename T>
std::vector<T> gttl_read_vector(const std::string &filename)
{
  return gttl_read_vector<T>(filename.c_str());
}
#endif
