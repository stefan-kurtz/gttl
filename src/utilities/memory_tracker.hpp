#ifndef MEMORY_TRACKER_HPP
#define MEMORY_TRACKER_HPP
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <map>
#include <cassert>
#include <cstdio>
#include "utilities/mathsupport.hpp"

#define GTTL_TRACK_MALLOC(TYPE,AMOUNT)\
        static_cast<TYPE *>(memory_tracker->malloc(__FILE__,__LINE__, AMOUNT))

#define GTTL_TRACK_CALLOC(TYPE,NUM_ELEMS,SIZE_ELEM)\
        static_cast<TYPE *>(memory_tracker->calloc(__FILE__,__LINE__,\
                                                   NUM_ELEMS,SIZE_ELEM))

#define GTTL_UNTRACK_ALLOC(PTR)\
        memory_tracker->untrack(PTR,__FILE__,__LINE__)

class GttlMemoryTracker
{
  struct Entry
  {
    const char *codefile_name;
    int codefile_line;
    size_t amount;
    Entry(void)
      : codefile_name(nullptr)
      , codefile_line(0)
      , amount(0)
    { }
    Entry(const char *_codefile_name,int _codefile_line,size_t _amount)
      : codefile_name(_codefile_name)
      , codefile_line(_codefile_line)
      , amount(_amount)
    { }
  };
  bool verbose;
  size_t current,
         max_so_far;
  std::map<void *,Entry> malloced_ptrs;
  public:
  GttlMemoryTracker(bool _verbose = false)
    : verbose(_verbose)
    , current(0)
    , max_so_far(0)
  { }
  ~GttlMemoryTracker(void)
  {
    for (auto &[ptr,entry] : malloced_ptrs)
    {
      fprintf(stderr,"memory leak of %zu bytes allocated in file %s in "
                     "line %d\n",
              entry.amount,
              entry.codefile_name,
              entry.codefile_line);
      exit(EXIT_FAILURE);
    }
    if (current > 0)
    {
      fprintf(stderr,"current = %zu > 0 not expected\n",current);
      exit(EXIT_FAILURE);
    }
  }
  void track(void *ptr,const char *codefile_name, int codefile_line,
             size_t amount)
  {
    current += amount;
    max_so_far = std::max(max_so_far,current);
    if (verbose)
    {
      printf("# memorytracker\t%s\t%d\t%zu\t%zu\n",
             codefile_name,codefile_line,current,max_so_far);
    }
    if (ptr != nullptr)
    {
      assert(not malloced_ptrs.contains(ptr));
      malloced_ptrs[ptr] = Entry(codefile_name,codefile_line,amount);
    }
  }
  void *malloc(const char *codefile_name, int codefile_line, size_t amount)
  {
    void *const ptr = std::malloc(amount);
    track(ptr,codefile_name,codefile_line,amount);
    return ptr;
  }
  void *calloc(const char *codefile_name, int codefile_line, size_t num_elems,
               size_t size_elem)
  {
    void *const ptr = std::calloc(num_elems, size_elem);
    track(ptr,codefile_name,codefile_line,num_elems * size_elem);
    return ptr;
  }
  void untrack(void *ptr,const char *codefile_name, int codefile_line)
  {
    if (not malloced_ptrs.contains(ptr))
    {
      fprintf(stderr,"file %s, line %d: cannot free memory\n",
                     codefile_name,codefile_line);
      exit(EXIT_FAILURE);
    }
    auto entry = malloced_ptrs[ptr];
    malloced_ptrs.erase(ptr);
    assert(current >= entry.amount);
    current -= entry.amount;
    if (verbose)
    {
      printf("# memoryuntrack\t%s\t%d\t%zu\t%zu\n",
             codefile_name,codefile_line,current,max_so_far);
    }
  }
  [[nodiscard]] size_t peak_get(void) const noexcept
  {
    return static_cast<size_t>(mega_bytes(max_so_far));
  }
  [[nodiscard]] size_t peak_in_bytes_get(void) const noexcept
  {
    return max_so_far;
  }
};
#endif
