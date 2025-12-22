#ifndef THREAD_POOL_VAR_HPP
#define THREAD_POOL_VAR_HPP
#include <cstdlib>
#include <cassert>
#include <thread>
#include <vector>
#include <memory>
#include "threading/virtual_queue.hpp"

class GttlThreadPoolVar
{
  public:
    template <class Fn, class... Args>
    GttlThreadPoolVar(size_t number_of_threads,
                      size_t number_of_tasks,
                      Fn && thread_func,
                      Args&&... args)
    {
      assert(number_of_threads >= 1 && number_of_tasks > 0);
      if (number_of_threads == 1)
      {
        for (size_t task_num = 0; task_num < number_of_tasks; task_num++)
        {
          thread_func(0, task_num, args...);
        }
      } else
      {
        std::vector<std::shared_ptr<std::thread>> threads{};
        VirtualQueue vq(number_of_tasks);
        threads.reserve(number_of_threads);
        for (size_t thd = 0; thd < number_of_threads; thd++)
        {
          threads.push_back(std::make_shared<std::thread>(
                             [&thread_func,&vq,thd,args...]() {
             size_t task_num;
             while ((task_num = vq.next_element()) <= vq.last_element())
             {
               thread_func(thd, task_num, args...);
             }
          }));
        }
        for (const auto& th : threads)
        {
          th->join();
        }
      }
    }
};
#endif
