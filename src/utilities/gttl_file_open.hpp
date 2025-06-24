#ifndef GTTL_FILE_OPEN_HPP
#define GTTL_FILE_OPEN_HPP

#include <cerrno>
#include <cstdio>
#include <string>
#include <stdexcept>
#include "utilities/file_size.hpp"
#include "utilities/has_gzip_header.hpp"

#ifndef GTTL_WITHOUT_ZLIB
#include <cstring>
#include <cassert>
#include <vector>
#include <zlib.h>

using GttlFpType = gzFile;
#define gttl_fp_type_open(FILENAME, MODE)   gzopen(FILENAME, MODE)
#define gttl_fp_type_close(FP)              gzclose(FP)
#define gttl_fp_type_reset(FP)              gzseek(FP, 0, SEEK_SET)
#define gttl_fp_type_gets(FP, BUFFER, SIZE) gzgets(FP, BUFFER, SIZE)
#define gttl_fp_type_getc(FP)               gzgetc(FP)
#define gttl_fp_type_is_eof(FP)             gzeof(FP)
#define gttl_fp_type_rewind(FP)             gzrewind(FP)
#define gttl_fp_type_read(BUF, SIZE, COUNT, FP) gzread(FP, BUF, (COUNT)*(SIZE))

#else
#include <ios>
using GttlFpType = FILE *;
#define gttl_fp_type_open(FILENAME, MODE)   std::fopen(FILENAME, MODE)
#define gttl_fp_type_close(FP)              fclose(FP)
#define gttl_fp_type_reset(FP)              fseek(FP, 0, SEEK_SET)
#define gttl_fp_type_gets(FP, BUFFER, SIZE) (fgets(BUFFER, SIZE, FP) != nullptr\
                                             ? BUFFER : nullptr)
#define gttl_fp_type_getc(FP)               getc(FP)
#define gttl_fp_type_is_eof(FP)             feof(FP)
#define gttl_fp_type_rewind(FP)             rewind(FP)
#define gttl_fp_type_read(BUF, SIZE, COUNT, FP) fread(BUF, SIZE, COUNT, FP)
#endif

template<typename BaseType>
static inline std::basic_string<BaseType>
  gttl_read_files(const std::vector<std::string> &inputfiles)
{
  std::vector<size_t> file_size_vec;
  file_size_vec.reserve(inputfiles.size());
  size_t sum_file_size = 0;
  bool append_sequences = false;
  for(auto &&inputfile : inputfiles)
  {
    if (has_gzip_header(inputfile.c_str()))
    {
#ifdef GTTL_WITHOUT_ZLIB
      throw std::ios_base::failure(std::string("cannot handle gzipped file ")
                                   + inputfile);
#else
      append_sequences = true;
      break;
#endif
    }
    const size_t file_size = gttl_file_size(inputfile.c_str());
    file_size_vec.push_back(file_size);
    sum_file_size += file_size;
  }
  std::basic_string<BaseType> concatenated_content;
  if (not append_sequences)
  {
    concatenated_content.resize(sum_file_size);
  }
  size_t offset = 0;
  size_t file_counter = 0;
  for (auto &&inputfile : inputfiles)
  {
    if (has_gzip_header(inputfile.c_str()))
    {
#ifndef GTTL_WITHOUT_ZLIB
      gzFile fp = gzopen(inputfile.c_str(), "rb");
      if (fp == nullptr)
      {
        throw std::runtime_error(std::string("Failed to gzopen file ")
                                 + inputfile);
      }
      assert(append_sequences);
      while (true)
      {
        constexpr const size_t buf_size = 4096;
        concatenated_content.resize(offset + buf_size);
        const int bytes_read = gzread(fp, concatenated_content.data() + offset,
                                      buf_size);
        if (bytes_read < 0)
        {
          throw std::runtime_error(std::string("Error reading from gzip file ")
                                   + inputfile);
        }
        assert (static_cast<size_t>(bytes_read) <= buf_size);
        if (static_cast<size_t>(bytes_read) < buf_size)
        {
          concatenated_content.resize(offset + bytes_read);
        }
        if (bytes_read == 0)
        {
          break;
        }
        offset += bytes_read;
      }
      gzclose(fp);
#endif
    } else
    {
      FILE *infp = std::fopen(inputfile.c_str(), "rb");
      if (infp == nullptr)
      {
        throw std::runtime_error(std::string("Error opening file: ")
                                 + std::to_string(errno));
      }
      const size_t file_size = append_sequences
                                 ? gttl_file_size(inputfile.c_str())
                                 : file_size_vec[file_counter++];
      if (append_sequences)
      {
        concatenated_content.resize(offset + file_size);
      }
      const size_t bytes_read = fread(concatenated_content.data() + offset,
                                      sizeof(BaseType), file_size, infp);
      if (bytes_read != file_size)
      {
        fclose(infp);
        throw std::runtime_error(std::string("Expected ")
                                 + std::to_string(file_size)
                                 + " bytes, but read only "
                                 + std::to_string(bytes_read));
      }
      fclose(infp);
      offset += file_size;
    }
  }
  return concatenated_content;
}

static inline std::string gttl_read_file(const char *file_name)
{
  const std::vector<std::string> inputfiles{std::string(file_name)};
  return gttl_read_files<char>(inputfiles);
}
#endif
