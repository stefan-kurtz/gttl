#ifndef OPTIONAL_LOCK_HPP
#define OPTIONAL_LOCK_HPP

#include <mutex>
class optional_lock
{
  std::mutex* mut;

  public:
  optional_lock(std::mutex* m) : mut(m)
  {
    if (mut != nullptr) mut->lock();
  }

  ~optional_lock()
  {
    if (mut != nullptr) mut->unlock();
  }
};

#endif // OPTIONAL_LOCK_HPP
