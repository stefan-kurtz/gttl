#ifndef GTTL_MMAP_HPP
#define GTTL_MMAP_HPP

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <ios>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include "utilities/file_size.hpp"
#include "utilities/str_format.hpp"

template<typename T>
class Gttlmmap
{
  int filedesc;
  size_t size_of_file,
         num_values;
  void *memorymap;
  public:
  Gttlmmap(const char *filename) :
    filedesc(0),
    size_of_file(gttl_file_size(filename)),
    num_values(size_of_file/sizeof(T)),
    memorymap(nullptr)
  {
    if (size_of_file % sizeof(T) != 0)
    {
      StrFormat msg("file %s contains %lu bytes which is not a multiple of %lu",
                    filename,size_of_file,sizeof(T));
      throw msg.str();
    }
    filedesc = open(filename, O_RDONLY);
    if (filedesc < 0)
    {
      StrFormat msg("cannot open file %s",filename);
      throw msg.str();
    }
    memorymap = mmap(NULL, size_of_file, PROT_READ, MAP_FILE | MAP_SHARED,
                     filedesc, 0);
    if (memorymap == MAP_FAILED)
    {
      StrFormat msg("cannot memory map %lu elements from file %s",
                    num_values,filename);
      throw msg.str();
    }
  }
  public:
  ~Gttlmmap(void)
  {
    assert(memorymap != nullptr);
    munmap(memorymap,size_of_file);
    assert(filedesc >= 0);
    close(filedesc);
  }
  const T *ptr(void) const noexcept
  {
    assert(memorymap != nullptr);
    return reinterpret_cast<const T *>(memorymap);
  }
  size_t size(void) const noexcept
  {
    return num_values;
  }
};
#endif
