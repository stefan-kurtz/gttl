#ifndef HAS_GZIP_HEADER_HPP
#define HAS_GZIP_HEADER_HPP
#include <cstdio>
#include <ios>
#include <string>
#include <cstring>

static inline bool has_gzip_header(const char *file_name)
{
  FILE *const fp = std::fopen(file_name, "rb");
  if (fp == nullptr)
  {
    throw std::ios_base::failure(": cannot open file "
                                 + std::string(file_name));
  }
  constexpr const unsigned char magic_bytes[] = {0x1F, 0x8B};
  unsigned char header[sizeof(magic_bytes)];
  if (fread(header, sizeof(unsigned char), sizeof(header), fp)
            != sizeof(header))
  {
    fclose(fp);
    // A file that does not have at least as many bytes as the GZip header
    // cannot be a valid GZip file by definition.
    return false;
  }
  fclose(fp);
  return memcmp(header, magic_bytes, sizeof(magic_bytes)) == 0;
}
#endif
