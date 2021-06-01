#ifndef ALL_VS_ALL_HPP
#define ALL_VS_ALL_HPP
#include <iostream>
#include <cstddef>
#include "threading/thread_pool.hpp"

template<typename T,typename R,typename Data,
         size_t (*first_index)(size_t),
         R (*process_pair)(const T &,size_t,const T &,size_t),
         void (*process_result)(size_t,size_t,size_t,R,Data &)>
static void gttl_one_vs_all(Data &data,size_t thread_id,
                            const T &tasks0,size_t idx,const T &tasks1)
{
  for (size_t j = first_index(idx); j < tasks1.size(); j++)
  {
    process_result(thread_id,idx,j,process_pair(tasks0,idx,tasks1,j),data);
  }
}

template<typename T,typename Data>
class GttlThreadData
{
  public:
    const T &tasks0, &tasks1;
    Data &data;
  GttlThreadData(const T &_tasks0,const T &_tasks1,Data &_d) :
    tasks0(_tasks0),
    tasks1(_tasks1),
    data(_d)
  {
  }
};

template<typename T,typename R,typename Data,size_t (*first_index)(size_t),
         R (*process_pair)(const T &,size_t,const T &,size_t),
         void (*process_result)(size_t,size_t,size_t,R,Data &)>
static void gttl_thread_runner(size_t thread_id,size_t task_num,
                               void *v_thread_data)
{
  GttlThreadData<T,Data> *thread_data
    = static_cast<GttlThreadData<T,Data> *>(v_thread_data);
  gttl_one_vs_all<T,R,Data,first_index,process_pair,process_result>
                      (thread_data->data,thread_id,
                       thread_data->tasks0,task_num,thread_data->tasks1);
}

template<typename T,typename R,typename Data,
         size_t (*first_index)(size_t),
         R (*process_pair)(const T &,size_t,const T &,size_t),
         void (*process_result)(size_t,size_t,size_t,R,Data &)>
void gttl_all_vs_all(const T &tasks0,const T &tasks1,Data &data,
                     size_t num_threads=1)
{
  assert(num_threads >= 1);
  if (num_threads == 1)
  {
    for (size_t idx = 0; idx < tasks0.size(); idx++)
    {
      gttl_one_vs_all<T,R,Data,first_index,process_pair,process_result>
                        (data,0,tasks0,idx,tasks1);
    }
  } else
  {
    GttlThreadData<T,Data> thread_data(tasks0,tasks1,data);
    GttlThreadPool thread_pool(num_threads,tasks0.size(),
                               gttl_thread_runner<T,R,Data,first_index,
                                                  process_pair,
                                                  process_result>,
                               (void *) &thread_data);
  }
}
#endif
