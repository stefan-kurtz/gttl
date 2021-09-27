#ifndef THREAD_POOL_HPP
#define THREAD_POOL_HPP
#include <cstdlib>
#include <cassert>
#include <thread>
#include <vector>
#include "virtual_queue.hpp"

template<class ThreadData>
using GttlThreadFunc = void (*)(size_t thread_id,size_t task_num,
                                ThreadData *thread_data);

template<class ThreadData>
static void gttl_thread_apply_thread_func(
                     GttlThreadFunc<ThreadData> thread_func,
                     ThreadData *thread_data,
                     VirtualQueue *vq,
                     size_t thread_id)
{
  size_t task_num;

  while ((task_num = vq->next_element()) <= vq->last_element())
  {
    thread_func(thread_id, task_num, thread_data);
  }
}

template<class ThreadData>
class GttlThreadPool
{
  public:
    GttlThreadPool(size_t number_of_threads,
                   size_t number_of_tasks,
                   GttlThreadFunc<ThreadData> thread_func,
                   ThreadData *thread_data)
    {
      assert(number_of_threads > 1 && number_of_tasks > 0);
      std::vector<std::thread> threads{};
      VirtualQueue vq(number_of_tasks);
      for (size_t thd = 0; thd < number_of_threads; thd++)
      {
        threads.push_back(std::thread(gttl_thread_apply_thread_func<ThreadData>,
                                      thread_func,
                                      thread_data,
                                      &vq, thd));
      }
      for (auto &th : threads)
      {
        th.join();
      }
    }
};
#endif
