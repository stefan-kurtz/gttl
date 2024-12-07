#ifndef GTTL_LINE_GENERATOR_COYIELD_HPP
#define GTTL_LINE_GENERATOR_COYIELD_HPP

#include <cstring>
#include <generator>
#include <string>
#include <vector>
#include "utilities/gttl_file_open.hpp"

template <const size_t buf_size = (1 << 14)>
std::generator<std::string_view> gttl_read_lines(GttlFpType fp)
{
  if (!fp)
  {
    throw std::ios_base::failure("Invalid file pointer");
  }

  // We use a vector instead of a static-sized array to prevent
  // problems with large stack allocation.
  static thread_local std::vector<char> buffer(buf_size);

  while (gttl_fp_type_gets(fp, buffer.data(), buf_size))
  {
    size_t len = std::strlen(buffer.data());
    // We trim LF and CR in the generator
    if (len > 0 && buffer[len - 1] == '\n')
    {
      buffer[len - 1] = '\0';
      len--;
      if (len > 0 && buffer[len - 1] == '\r')
      {
        buffer[len - 1] = '\0';
        len--;
      }
    }
    //We return a string_view directly into the buffer
    co_yield std::string_view(buffer.data(), len);
  }

  gttl_fp_type_close(fp);
}

template <const size_t buf_size = (1 << 14)>
std::generator<std::string> gttl_read_lines(const std::string file_path)
{
  return gttl_read_lines<buf_size>(gttl_fp_type_open(file_path.c_str(), "r"));
}

#endif  // GTTL_LINE_GENERATOR_COYIELD_HPP
