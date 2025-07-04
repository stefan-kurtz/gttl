#ifndef WRITE_OUTPUT_FILE_HPP
#define WRITE_OUTPUT_FILE_HPP

#include <zlib.h>
#include <cerrno>
#include <cstdio>
#include <ios>
#include <string>
#include <system_error>

/*
** A helper function to write (GZip) compressed or uncompressed output-files.
** compression_level=0 will result in uncompressed output.
** .gz file extensions will be automatically added where appropriate.
 */
static inline void write_to_output_file(const std::string &file_name,
                                 const std::string &content,
                                 size_t compression_level)
{
  if(compression_level > 9)
  {
    throw std::ios_base::failure(
      ": GZip compression level can not be greater than 9!");
  }
  if(compression_level == 0)
  {
    FILE *const f_out = std::fopen(file_name.c_str(), "w");
    if(f_out == nullptr)
    {
      throw std::ios_base::failure(": Error writing to file: " + file_name);
    }
    fputs(content.c_str(), f_out);
    fclose(f_out);
  } else
  {
    const std::string gz_file_name = file_name + ".gz";
    const std::string write_mode = "w" + std::to_string(compression_level);
    const gzFile f_out = gzopen(gz_file_name.c_str(), write_mode.c_str());
    if(f_out == Z_NULL)
    {
      throw std::system_error(errno,
                              std::iostream_category(),
                              ": Error writing to file: " + gz_file_name);
    }
    gzputs(f_out, content.c_str());
    gzclose(f_out);
  }
}

#endif // WRITE_OUTPUT_FILE_HPP
