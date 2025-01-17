#ifndef GTTL_FILE_OPEN_HPP
#define GTTL_FILE_OPEN_HPP

#ifndef GTTL_WITHOUT_ZLIB
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <zlib.h>
#include "utilities/file_size.hpp"
using GttlFpType = gzFile;
#define gttl_fp_type_open(FILENAME, MODE)   gzopen(FILENAME, MODE)
#define gttl_fp_type_close(FP)              gzclose(FP)
#define gttl_fp_type_reset(FP)              gzseek(FP, 0, SEEK_SET)
#define gttl_fp_type_gets(FP, BUFFER, SIZE) gzgets(FP, BUFFER, SIZE)
#define gttl_fp_type_is_eof(FP)             gzeof(FP)

static inline std::string gttl_read_file(const char *file_name)
{
  FILE* fp = fopen(file_name, "rb");
  if(!fp)
  {
    throw std::runtime_error("Failed to open file: " + std::string(file_name));
  }

  constexpr const unsigned char magic_bytes[] = {0x1F, 0x8B};
  unsigned char header[sizeof(magic_bytes)];

  if(fread(header, sizeof(unsigned char), sizeof(header), fp) != sizeof(header))
  {
    fclose(fp);
    throw std::runtime_error("Failed to read header of file: "
                             + std::string(file_name));
  }

  // Check if the file is compressed by comparing magic-bytes
  if(memcmp(header, magic_bytes, sizeof(magic_bytes)) != 0)
  {
    // Not compressed
    rewind(fp);
    size_t file_size = gttl_file_size(file_name);

    std::string content(file_size, '\0');
    if(fread(content.data(), 1, file_size, fp) != file_size)
    {
      fclose(fp);
      throw std::runtime_error("Failed to read entire file: "
                               + std::string(file_name));
    }

    fclose(fp);
    return content;
  }else
  {
    fclose(fp);
    gzFile fp = gzopen(file_name, "rb");
    if(!fp)
    {
      throw std::runtime_error("Failed to gzopen file: "
                               + std::string(file_name));
    }
    constexpr const size_t buf_size = 4096;
    std::string content;
    int bytes_read;
    char buffer[buf_size];

    while((bytes_read = gzread(fp, buffer, buf_size)))
    {
      content.append(buffer, bytes_read);
    }

    gzclose(fp);

    if(bytes_read < 0)
    {
      throw std::runtime_error(std::string("Error reading from gzip file: ")
                               + file_name);
    }

    return content;
  }
}

#else
using GttlFpType = FILE *;
#include <cstdio>
#include "utilities/file_size.hpp"
#define gttl_fp_type_open(FILENAME, MODE)   fopen(FILENAME, MODE)
#define gttl_fp_type_close(FP)              fclose(FP)
#define gttl_fp_type_reset(FP)              fseek(FP, 0, SEEK_SET)
#define gttl_fp_type_gets(FP, BUFFER, SIZE) (fgets(BUFFER, SIZE, FP) != nullptr\
                                             ? BUFFER : nullptr)
#define gttl_fp_type_is_eof(FP)             feof(FP)

static inline std::string gttl_read_file(const char *file_name)
{
  std::string content;
  content.resize(gttl_file_size(file_name));
  FILE* fp = fopen(file_name, "rb");
  if(fp)
  {
    fread(content.data(), sizeof(char), content.size(), fp);
    return content;
  }
  throw(errno);
}
#endif
#endif
