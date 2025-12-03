#ifndef THREADSAFE_QUEUE_HPP
#define THREADSAFE_QUEUE_HPP
#include <algorithm>
#include <exception>
#include <cstddef>
#include <mutex>
#include <optional>
#include <queue>
#include <string>
#include <utility>

// Source from https://bitbucket.org/marco/samples.git queue.cpp

class non_empty_queue : public std::exception
{
  std::string what_;
 public:
  explicit non_empty_queue(std::string msg) : what_(std::move(msg)) {}
  [[nodiscard]] const char* what(void) const noexcept override
  {
    return what_.c_str();
  }
};

template<typename T>
class ThreadsafeQueue
{
  std::queue<T> the_queue;
  mutable std::mutex q_mutex;

  // not public to prevent races between this and dequeue().
  bool empty() const
  {
    return the_queue.empty();
  }

 public:
  ThreadsafeQueue(void) = default;
  ThreadsafeQueue(const ThreadsafeQueue<T> &) = delete ;
  ThreadsafeQueue& operator=(const ThreadsafeQueue<T> &) = delete ;

  ThreadsafeQueue(ThreadsafeQueue<T>&& other) noexcept(false)
  {
    std::scoped_lock<std::mutex> lock(q_mutex);
    if (not empty())
    {
      throw non_empty_queue(std::string("Moving into a non-empty queue"));
    }
    the_queue = std::move(other.the_queue);
  }

  ~ThreadsafeQueue(void) noexcept(false)
  {
    const std::scoped_lock<std::mutex> lock(q_mutex);
    if (not empty())
    {
      throw non_empty_queue(std::string("Destroying a non-empty queue"));
    }
  }


  void enqueue(const T &item)
  {
    const std::scoped_lock<std::mutex> lock(q_mutex);
    the_queue.push(item);
  }

  std::optional<T> dequeue(void)
  {
    const std::scoped_lock<std::mutex> lock(q_mutex);
    if (the_queue.empty())
    {
      return {};
    }
    T tmp = the_queue.front();
    the_queue.pop();
    return tmp;
  }

  size_t size(void) const
  {
    const std::scoped_lock<std::mutex> lock(q_mutex);
    return the_queue.size();
  }

};
#endif
