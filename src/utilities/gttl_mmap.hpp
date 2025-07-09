#ifndef GTTL_MMAP_HPP
#define GTTL_MMAP_HPP

#include <iostream>
#include <cstdlib>
#include <ios>
#include <cassert>
#include <fcntl.h>
#include "utilities/file_size.hpp"
#include "utilities/str_format.hpp"
#ifdef _WIN32
  #define NOMINMAX
  #include <io.h>
  #include "utilities/windows_mman.hpp"
#else
  #include <unistd.h>
  #include <sys/mman.h>
#endif

template<typename T>
class Gttlmmap
{
  int filedesc;
  size_t size_of_file,
         num_values;
  void *memorymap;
  public:
  explicit Gttlmmap(const char *filename) :
    filedesc(0),
    size_of_file(gttl_file_size(filename)),
    num_values(size_of_file/sizeof(T)),
    memorymap(nullptr)
  {
    if (size_of_file % sizeof(T) != 0)
    {
      const StrFormat msg("file %s contains %zu bytes which is not a "
                          "multiple of %zu",
                          filename,
                          size_of_file,
                          sizeof(T));
      throw std::ios_base::failure(msg.str());
    }
    filedesc = open(filename, O_RDONLY);
    if (filedesc < 0)
    {
      const StrFormat msg("cannot open file %s", filename);
      throw std::ios_base::failure(msg.str());
    }
    memorymap = mmap(nullptr, size_of_file, PROT_READ, MAP_FILE | MAP_SHARED,
                     filedesc, 0);
    if (memorymap == MAP_FAILED)
    {
      const StrFormat msg("cannot memory map %zu elements from file %s",
                          num_values,
                          filename);
      throw std::ios_base::failure(msg.str());
    }
  }

  ~Gttlmmap(void)
  {
    assert(memorymap != nullptr);
    munmap(memorymap,size_of_file);
    assert(filedesc >= 0);
    close(filedesc);
  }
  [[nodiscard]] const T *ptr(void) const noexcept
  {
    assert(memorymap != nullptr);
    return reinterpret_cast<const T *>(memorymap);
  }
  [[nodiscard]] size_t size(void) const noexcept
  {
    return num_values;
  }
};
#endif
