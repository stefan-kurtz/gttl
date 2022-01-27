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

class GttlThreadPool
{
  public:
    template<class ThreadData>
    GttlThreadPool(size_t number_of_threads,
                   size_t number_of_tasks,
                   GttlThreadFunc<ThreadData> thread_func,
                   ThreadData *thread_data)
    {
      assert(number_of_threads >= 1 && number_of_tasks > 0);
      if (number_of_threads == 1)
      {
        for (size_t task_num = 0; task_num < number_of_tasks; task_num++)
        {
          thread_func(0, task_num, thread_data);
        }
      } else
      {
        std::vector<std::thread> threads{};
        VirtualQueue vq(number_of_tasks);
        for (size_t thd = 0; thd < number_of_threads; thd++)
        {
          threads.push_back(std::thread([&thread_func,&thread_data,&vq,thd]() {
             size_t task_num;
             while ((task_num = vq.next_element()) <= vq.last_element())
             {
               thread_func(thd, task_num, thread_data);
             }
          }));
        }
        for (auto &th : threads)
        {
          th.join();
        }
      }
    }
};
#endif
