#ifndef OPTIONAL_LOCK_HPP
#define OPTIONAL_LOCK_HPP

#include <mutex>
class OptionalLock
{
  std::mutex* mut;

  public:
  OptionalLock(std::mutex* m) : mut(m)
  {
    if (mut != nullptr) mut->lock();
  }

  ~OptionalLock()
  {
    if (mut != nullptr) mut->unlock();
  }
};

#endif // OPTIONAL_LOCK_HPP
