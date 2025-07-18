#ifndef TRACKING_ALLOCATOR_HPP
#define TRACKING_ALLOCATOR_HPP

#include <atomic>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <type_traits>
#include "utilities/mathsupport.hpp"

#if defined(__clang__) || defined(__GNUC__)

#include <cxxabi.h>
#include <string>

inline std::string demangle(const char* name)
{
  // Written by Ali Baharev in
  // https://stackoverflow.com/questions/281818/
  // Comments added by Stevan Miladinovic
  int status = 0;
  // abi::__cva_demangle returns a char* to memory
  // allocated with std::malloc().
  // This means that we have to explicitly free it ourselves, using
  // std::free().
  // A std::unique_ptr takes care of automatically freeing memory as the
  // last reference is lost, but it defaults to delete, not std::free().
  // Hence why we need to pass std::free as a "custom" deallocator here.
  const std::unique_ptr<char, void (*)(void *)> demangled(
                               abi::__cxa_demangle(name,
                                                   nullptr,
                                                   nullptr,
                                                   &status),
                               std::free);
  return (status == 0 and demangled) ? demangled.get() : name;
}

#else

// On MSVC, these should be readable. Elsewhere, this is yet to be implemented.
inline std::string demangle(const char* name)
{
  return name;
}

#endif

template <typename T>
struct TrackingAllocator : std::allocator<T>
{
  using value_type = T;
  static_assert(not std::is_void_v<T>,
                "TrackingAllocator does not support void-types!\n");

  static inline std::atomic_size_t total_allocated = 0;
  static inline std::atomic_size_t max_allocated = 0;
  static inline std::atomic_size_t allocations = 0;
  static inline std::atomic_size_t deallocations = 0;

  T* allocate(size_t n)
  {
    // Using a separate variable may reduce atomic-overhead
    // on low optimization levels
    const size_t bytes = n * sizeof(T);
    total_allocated += bytes;

    const size_t current = total_allocated.load();
    size_t prev_max = max_allocated.load();
    // This loop simply updates the max-value in a thread-safe way.
    // The variable itself is atomic, but something might change between
    // the comparison (atomic read) and the update (atomic write).
    // Hence we need to explicitly re-check with compare_exchange_weak.
    while (current > prev_max
           and not max_allocated.compare_exchange_weak(prev_max, current))
    {}
    ++allocations;
    return std::allocator<T>::allocate(n);
  }

  void deallocate(T* p, size_t n)
  {
    total_allocated -= n * sizeof(T);
    ++deallocations;
    std::allocator<T>::deallocate(p, n);
  }

  // This struct is used when a container allocated using
  // TrackingAllocator<T> needs to allocate memory itself.
  // A good example of this is a std::unordered_set<T>.
  // Here, a bucket<T> is first allocated in the set, into which Ts are then
  // placed.
  // This would then call TrackingAllocator<T>::rebind<bucket<T>> to allocate
  // the bucket.
  template <typename U>
  struct rebind
  {
    using other = TrackingAllocator<U>;
  };

  // Since all members are static, there is effectively
  // only one tracking allocator
  bool operator==(const TrackingAllocator& /**/) const noexcept
  { return true; }
  bool operator!=(const TrackingAllocator& /**/) const noexcept
  { return false; }

  static void log_stats()
  {
    printf("# TrackingAllocator<%s> "
           "max_allocated\t%.1f\t"
           "allocations\t%zu\t"
           "deallocations\t%zu\n",demangle(typeid(T).name()).c_str(),
                                  mega_bytes(max_allocated.load()),
                                  allocations.load(),
                                  deallocations.load());
    if(allocations != deallocations)
    {
      fprintf(stderr,"WARNING: #allocations != #deallocations "
                     "=> memory leak\n");
    }
  }

  TrackingAllocator() noexcept
  {
    register_atexit_logger();
  };

  // We need this to allow rebinding via the above struct.
  // to continue the example from above, a TrackingAllocator<bucket<T>> needs
  // the be able to then construct as many TrackingAllocator<T> as required.
  // It will pass itself as a reference in this process. Though we do not
  // need this argument for our use-case, we still need to accept
  // it in a "copy-constructor" of sorts to satisfy STL requirements.
  template <typename U>
  explicit TrackingAllocator(const TrackingAllocator<U>& /**/) noexcept
  {
    register_atexit_logger();
  }

  private:
  static void register_atexit_logger()
  {
    static const bool registered = ([]
    {
      // std::atexit simply runs as the program exits, including on potential
      // SIGTERM.
      if(std::atexit(log_stats) != 0)
      {
        fprintf(stderr,"WARNING: #allocations != #deallocations "
                       "=> memory leak\n");
        // We don't throw, because this is not critical
      }
      return true;
    })();
    (void)registered;
  }
};


#endif // TRACKING_ALLOCATOR_HPP
