#ifndef FILE_SIZE_HPP
#define FILE_SIZE_HPP

/* Use some C based functions to determine size of file in portable
   way, as g++ 7.5 and clang++ on Linux do not know <filesystem>. */

#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>  // for fstat() and open()
#include <sys/stat.h>
#include "utilities/str_format.hpp"

inline size_t gttl_file_size(const char *filename)
{
  int filedesc = open(filename,O_RDONLY);
  if (filedesc == -1)       // check for error code
  {
    StrFormat msg("cannot open file %s",filename);
    throw msg.str();
  }
  struct stat buf;
  if (fstat(filedesc,&buf) == -1) // get status of file
  {
    StrFormat msg("cannot access status of file %s",filename);
    throw msg.str();
  }
  close(filedesc);
  return static_cast<size_t>(buf.st_size);
}
#endif
