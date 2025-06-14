#ifndef TRACKING_ALLOCATOR_HPP
#define TRACKING_ALLOCATOR_HPP

#include <atomic>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <type_traits>

#if defined(__clang__) || defined(__GNUC__)

#include <cxxabi.h>
#include <string>

inline std::string demangle(const char* name)
{
  int status = 0;
  std::unique_ptr<char, void(*)(void*)> demangled(
    abi::__cxa_demangle(name, nullptr, nullptr, &status),
    std::free
  );
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
    std::cerr << "TrackingAllocator<" << demangle(typeid(T).name())
              << "> stats:\n"
              << "  max_allocated:\t" << max_allocated.load() << " bytes\n"
              << "  allocations:\t"   << allocations.load() << "\n"
              << "  deallocations:\t" << deallocations.load() << "\n";

    if(allocations != deallocations)
    {
      std::cerr << "WARNING: allocations != deallocations!\n"
                << "This indicates a memory leak!\n";
    }
  }

  TrackingAllocator() noexcept
  {
    register_atexit_logger();
  };

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
      if(std::atexit(log_stats) != 0)
      {
        std::cerr << "Warning: Failed to register the"
                  << "TrackingAllocator logger with atexit()!\n";
        // We don't throw, because this is not critical
      }
      return true;
    })();
    (void)registered;
  }
};


#endif // TRACKING_ALLOCATOR_HPP
