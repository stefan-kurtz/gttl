#ifndef THREAD_POOL_UNKNOWN_TASKS_HPP
#define THREAD_POOL_UNKNOWN_TASKS_HPP

#include <condition_variable>
#include <mutex>
#include <queue>
#include <thread>
#include "threading/threadsafe_queue.hpp"


template <class FunctionType>
class ThreadPoolUnknownTasks
{
  private:
    std::vector<std::thread> threads;
    std::queue<FunctionType> tasks;
    ThreadsafeQueue<FunctionType> tsq;
    std::mutex queue_lock;
    std::condition_variable tasks_changed;
    bool stop = false;
  public:
    ThreadPoolUnknownTasks(const size_t num_threads =
                           std::thread::hardware_concurrency())
    {
      for(size_t i = 0; i < num_threads; i++)
      {
        threads.emplace_back([this]
        {
          while(true)
          {
            FunctionType task;
            {
              std::unique_lock<std::mutex> lock(queue_lock);
              tasks_changed.wait(lock, [this]{return (tsq.size()!=0) || stop;});
              if(stop && tsq.size() == 0)
              {
                return;
              }
              task = *(tsq.dequeue());
            }
            task();
          }
        });
      }
    }

    ~ThreadPoolUnknownTasks()
    {
      {
        std::unique_lock<std::mutex> lock(queue_lock);
        stop = true;
      }
      tasks_changed.notify_all();
      for(auto& t : threads)
      {
        t.join();
      }
    }

    void enqueue(FunctionType task)
    {
      tsq.enqueue(task);
      tasks_changed.notify_one();
    }
};


#endif // THREAD_POOL_UNKNOWN_TASKS_HPP
