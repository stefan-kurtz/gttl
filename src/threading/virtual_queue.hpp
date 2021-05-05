#ifndef ATOMIC_QUEUE_HPP
#define ATOMIC_QUEUE_HPP
#include <cstdlib>
#include <atomic>

class VirtualQueue {
  private:
    std::atomic<size_t> current;
    size_t last;
  public:
    VirtualQueue(size_t number_of_elements) {
      current = 0;
      assert(number_of_elements > 0);
      last = number_of_elements - 1;
    }
    ~VirtualQueue(void) {}
    size_t next_element(void) {
      return current++;
    }
    size_t last_element(void) const {
      return last;
   }
};
#endif
