#ifndef RUNTIME_CLASS_HPP
#define RUNTIME_CLASS_HPP
#include <iostream>
#include <mutex>
#include <cassert>
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
  std::string to_string(void)
  {
    auto end_time = std::chrono::high_resolution_clock::now();
    size_t elapsed_micro
      = (size_t) std::chrono::duration_cast<std::chrono::microseconds>
                      (end_time-start_time).count();
    std::string s = std::to_string(this->to_ms(elapsed_micro));
    start_time = end_time;
    return s;
  }
  std::string to_string(const char *msg)
  {
    return std::string("TIME\t") + std::string(msg) + std::string(" (ms):\t")
                                 + this->to_string();
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
    assert(cout_mutex != nullptr);
    cout_mutex->lock();
    const size_t elapsed_micro = this->show(msg);
    cout_mutex->unlock();
    return elapsed_micro;
  }
  size_t locked_show(std::mutex *cout_mutex,const std::string &msg)
  {
    assert(cout_mutex != nullptr);
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
