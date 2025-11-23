#ifndef VIRTUAL_QUEUE_HPP
#define VIRTUAL_QUEUE_HPP
#include <cassert>
#include <cstdlib>
#include <atomic>

class VirtualQueue {
  private:
    std::atomic<size_t> current;
    size_t last;
  public:
    VirtualQueue(size_t number_of_elements)
      : current(0)
      , last(number_of_elements - 1)
    {
      assert(number_of_elements > 0);
    }
    ~VirtualQueue(void) = default;
    size_t next_element(void)
    {
      return current++;
    }
    [[nodiscard]] size_t last_element(void) const noexcept { return last; }
};
#endif
