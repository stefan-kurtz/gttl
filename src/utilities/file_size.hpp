#ifndef FILE_SIZE_HPP
#define FILE_SIZE_HPP

/* Use some C based functions to determine size of file in portable
   way, as g++ 7.5 and clang++ on Linux do not know <filesystem>. */

#ifdef _WIN32
  #define NOMINMAX
  #include <io.h>
#else
  #include <unistd.h>
#endif
#include <cstddef>
#include <fcntl.h>
#include <ios>
#include <sys/stat.h>
#include <string>
#include <vector>

inline size_t gttl_file_size(const char *filename)
{
  const int filedesc = open(filename, O_RDONLY);
  if (filedesc == -1)       // check for error code
  {
    throw std::ios_base::failure(": cannot open file " + std::string(filename));
  }
  struct stat buf;
  if (fstat(filedesc,&buf) == -1) // get status of file
  {
    throw std::ios_base::failure(
          std::string(": cannot access status of file ") +
          std::string(filename));
  }
  close(filedesc);
  return static_cast<size_t>(buf.st_size);
}

inline size_t gttl_file_size(const std::string &inputfile)
{
  return gttl_file_size(inputfile.c_str());
}

inline size_t gttl_file_size(const std::vector<std::string> &inputfiles)
{
  size_t files_bytes = 0;
  for (const auto &inputfile : inputfiles)
  {
    files_bytes += gttl_file_size(inputfile);
  }
  return files_bytes;
}
#endif
