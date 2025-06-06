#ifndef GTTL_BINARY_WRITE_HPP
#define GTTL_BINARY_WRITE_HPP
#include <cstddef>
#include <cstdio>
#include <cassert>
#include <ios>
#include <string>

template <typename T,
          size_t buf_size = (size_t(1) << 16)/ sizeof(T)>
class BinaryFileWriter
{
  static_assert(buf_size > 0, "Buffer size may not be zero");
  static_assert(std::is_trivially_copyable_v<T>,
                "BinaryFileWriter can only work with types that are "
                "trivially copyable.");
  T buffer[buf_size];
  FILE *out_fp;
  size_t nextfree = 0;
  public:
  BinaryFileWriter(const std::string &outfilename)
    : out_fp(fopen(outfilename.c_str(), "wb"))
  {
    if (out_fp == nullptr)
    {
      throw std::ios_base::failure(
              std::string(": cannot create file \"") +
              std::string(outfilename) +
              std::string("\""));
    }
  }
  void append(T value)
  {
    if (nextfree == buf_size)
    {
      std::fwrite(buffer, sizeof(T), buf_size,  out_fp);
      nextfree = 0;
    }
    assert(nextfree < buf_size);
    buffer[nextfree++] = value;
  }
  ~BinaryFileWriter(void)
  {
    if (nextfree > 0)
    {
      std::fwrite(buffer, sizeof(T), nextfree,  out_fp);
    }
    fclose(out_fp);
  }
};
#endif
