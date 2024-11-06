#ifndef WRITE_OUTPUT_FILE_HPP
#define WRITE_OUTPUT_FILE_HPP

#include <zlib.h>
#include <cstdio>
#include <string>
#include "utilities/str_format.hpp"

/*
** A helper function to write (GZip) compressed or uncompressed output-files.
** compression_level=0 will result in uncompressed output.
** .gz file extensions will be automatically added where appropriate.
 */
static void write_to_output_file(std::string file_name,
                                 const std::string &content,
                                 size_t compression_level)
{
  if(compression_level > 9)
  {
    StrFormat msg(": GZip compression level can not be greater than 9!");
    throw msg.str();
  }
  if(compression_level == 0)
  {
    FILE *f_out = fopen(file_name.c_str(), "w");
    if(f_out == NULL)
    {
      StrFormat msg(": Error writing to file: %s", file_name.c_str());
      throw msg.str();
    }

    fputs(content.c_str(), f_out);
    fclose(f_out);
  }else
  {
    file_name += ".gz";

    const std::string write_mode = "w"+std::to_string(compression_level);
    gzFile f_out = gzopen(file_name.c_str(), write_mode.c_str());
    if(f_out == Z_NULL)
    {
      StrFormat msg(": Error writing to file: %s", file_name.c_str());
      throw msg.str();
    }

    gzputs(f_out, content.c_str());
    gzclose(f_out);
  }
}

#endif // WRITE_OUTPUT_FILE_HPP
