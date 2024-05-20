#ifndef SAFE_QUEUE
#define SAFE_QUEUE

#include <condition_variable>
#include <mutex>
#include <queue>

// A threadsafe-queue from
// https://gist.githubusercontent.com/murphypei/safe_queue.h
template <class T>
class SafeQueue
{
  private:
    std::queue<T> this_queue;
    mutable std::mutex this_mutex;
    std::condition_variable cv;
  public:
    SafeQueue()
      : this_queue()
      , this_mutex()
      , cv()
    {}

    // Add an element to the queue.
    void enqueue(T t)
    {
      std::lock_guard<std::mutex> lock(this_mutex);
      this_queue.push(t);
      cv.notify_one();
    }

    // Get the front element.
    // If the queue is empty, wait till a element is avaiable.
    T dequeue(void)
    {
      std::unique_lock<std::mutex> lock(this_mutex);
      while (this_queue.empty())
      {
        // release lock as long as the wait and reaquire it afterwards.
        cv.wait(lock);
      }
      T val = this_queue.front();
      this_queue.pop();
      return val;
    }
};
#endif
