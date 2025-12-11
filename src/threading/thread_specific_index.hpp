#ifndef THREAD_SPECIFIC_INDEX_HPP
#define THREAD_SPECIFIC_INDEX_HPP
#include <cstddef>
#include <format>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <map>

class ThreadSpecificIndex
{
  const size_t num_threads;
  std::mutex thread_id_mutex{};
  std::map<std::thread::id, size_t> thread_id_map;
  public:
  ThreadSpecificIndex(size_t _num_threads)
    : num_threads(_num_threads)
  { }
  size_t get(void)
  {
    const std::thread::id t_id = std::this_thread::get_id();
    const std::scoped_lock<std::mutex> thread_id_lock(thread_id_mutex);
    if (thread_id_map.size() < num_threads)
    {
      const size_t current_size = thread_id_map.size();
      thread_id_map[t_id] = current_size;
    } else
    {
      if (thread_id_map.size() != num_threads)
      {
        throw std::runtime_error(std::format("thread_id_map.size() = {} != "
                                             "{} = num_threads",
                                             thread_id_map.size(),
                                             num_threads));
      }
    }
    return thread_id_map[t_id];
  }
};
#endif
