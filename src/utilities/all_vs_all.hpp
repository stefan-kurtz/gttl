#include <iostream>
#include <cstddef>
#include "threading/thread_pool.hpp"

template<typename T,typename R,typename Data,
         R (*process_pair)(T &,T &),
         void (*process_result)(size_t,size_t,size_t,R,Data &)>
static void gttl_one_vs_all(Data &data,size_t thread_id,
                            size_t idx,std::vector<T> &tasks)
{
  for (size_t j = idx+1; j < tasks.size(); j++)
  {
    process_result(thread_id,idx,j,process_pair(tasks[idx],tasks[j]),data);
  }
}

template<typename T,typename Data>
class GttlThreadData
{
  public:
    std::vector<T> &tasks;
    Data &data;
  GttlThreadData(std::vector<T> &_t,Data &_d) :
    tasks(_t),
    data(_d)
  {
  }
};

template<typename T,typename R,typename Data,R (*process_pair)(T &,T &),
         void (*process_result)(size_t,size_t,size_t,R,Data &)>
static void gttl_thread_runner(size_t thread_id,size_t task_num,
                               void *v_thread_data)
{
  GttlThreadData<T,Data> *thread_data
    = static_cast<GttlThreadData<T,Data> *>(v_thread_data);
  gttl_one_vs_all<T,R,Data,process_pair,process_result>(thread_data->data,
                                                        thread_id,
                                                        task_num,
                                                        thread_data->tasks);
}

template<typename T,typename R,typename Data,R (*process_pair)(T &,T &),
         void (*process_result)(size_t,size_t,size_t,R,Data &)>
void gttl_all_vs_all(std::vector<T> &tasks,Data &data,size_t num_threads=1)
{
  assert(num_threads >= 1);
  if (num_threads == 1)
  {
    for (size_t idx = 0; idx < tasks.size(); idx++)
    {
      gttl_one_vs_all<T,R,Data,process_pair,process_result>(data,0,idx,tasks);
    }
  } else
  {
    GttlThreadData<T,Data> thread_data(tasks,data);
    GttlThreadPool thread_pool(num_threads,tasks.size(),
                               gttl_thread_runner<T,R,Data,process_pair,
                                                  process_result>,
                               (void *) &thread_data);
  }
}
