#ifndef RUN_TIME_CLASS_HPP
#define RUN_TIME_CLASS_HPP
#include <iostream>
#include <mutex>
#include <chrono>

class RunTimeClass
{
  private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
  public:
  RunTimeClass(void)
  {
    start_time = std::chrono::high_resolution_clock::now();
  }
  size_t to_ms(size_t time_microseconds)
  {
    return time_microseconds/size_t(1000);
  }
  size_t show(const char *msg)
  {
    auto end_time = std::chrono::high_resolution_clock::now();
    size_t elapsed_micro
      = (size_t) std::chrono::duration_cast<std::chrono::microseconds>
                      (end_time-start_time).count();
    std::cout << "# TIME\t" << msg << " (ms):\t" << this->to_ms(elapsed_micro)
              << std::endl;
    start_time = end_time;
    return elapsed_micro;
  }
  size_t show(const std::string &msg)
  {
   return this->show(msg.c_str());
  }
  size_t locked_show(std::mutex *cout_mutex,const char *msg)
  {
    cout_mutex->lock();
    const size_t elapsed_micro = this->show(msg);
    cout_mutex->unlock();
    return elapsed_micro;
  }
  size_t elapsed(void)
  {
    auto end_time = std::chrono::high_resolution_clock::now();
    return (size_t) std::chrono::duration_cast<std::chrono::microseconds>
                         (end_time-start_time).count();
  }
  void reset(void)
  {
    start_time = std::chrono::high_resolution_clock::now();
  }
};
#endif
