#ifndef FILE_EXISTS_HPP
#define FILE_EXISTS_HPP
#include <sys/stat.h>
#include <cstdbool>
#include <string>
static inline bool gttl_file_exists(const std::string &filename)
{
  struct stat buf;
  return stat (filename.c_str(), &buf) == 0;
}
#endif
