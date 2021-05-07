#ifndef THREAD_POOL_HPP
#define THREAD_POOL_HPP
#include <cstdlib>
#include <cassert>
#include <thread>
#include "virtual_queue.hpp"

typedef void (*GttlThreadFunc)(size_t thread_id,size_t task_num,
                               void *thread_data);

static void gttl_thread_apply_thread_func(GttlThreadFunc thread_func,
                                          void *thread_data,
                                          VirtualQueue *vq,
                                          size_t thread_id)
{
  size_t task_num;

  while ((task_num = vq->next_element()) <= vq->last_element())
  {
    thread_func(thread_id, task_num, thread_data);
  }
}

class GttlThreadPool
{
  public:
  GttlThreadPool(size_t num_threads, size_t number_of_tasks,
                 GttlThreadFunc thread_func,void *thread_data)
  {
    assert(num_threads >= 1 && number_of_tasks > 0);
    if (num_threads > 1)
    {
      std::thread *threads = new std::thread [num_threads];
      VirtualQueue vq(number_of_tasks);
      for (size_t thd = 0; thd < num_threads; thd++)
      {
        threads[thd] = std::thread(gttl_thread_apply_thread_func,
                                   thread_func,thread_data,
                                   &vq, thd);
      }
      for (size_t thd = 0; thd < num_threads; thd++)
      {
        threads[thd].join();
      }
      delete[] threads;
    } else
    {
      for (size_t task_num = 0; task_num < number_of_tasks; task_num++)
      {
        thread_func(0,task_num,thread_data);
      }
    }
  }
};
#endif
