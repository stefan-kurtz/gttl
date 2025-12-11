#ifndef RUNTIME_CLASS_HPP
#define RUNTIME_CLASS_HPP
#include <cstdio>
#include <mutex>
#include <cassert>
#include <chrono>
#include <string>

class RunTimeClass
{
  private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
  public:
  RunTimeClass(void)
  {
    start_time = std::chrono::high_resolution_clock::now();
  }
  [[nodiscard]] size_t to_ms(size_t time_microseconds) const
  {
    return time_microseconds/size_t(1000);
  }
  std::string to_string(void)
  {
    const auto end_time = std::chrono::high_resolution_clock::now();
    const size_t elapsed_micro
      = static_cast<size_t>
                   (std::chrono::duration_cast<std::chrono::microseconds>
                    (end_time - start_time).count());
    std::string s = std::to_string(this->to_ms(elapsed_micro));
    start_time = end_time;
    return s;
  }
  std::string to_string(const std::string &msg)
  {
    return std::string("TIME\t") + msg + std::string(" (ms):\t")
                                 + this->to_string();
  }
  std::string to_string(const char *msg)
  {
    return to_string(std::string(msg));
  }
  size_t show_fp(FILE *out_fp, const char *msg)
  {
    auto end_time = std::chrono::high_resolution_clock::now();
    const size_t elapsed_micro = (size_t) std::chrono::duration_cast<
                                                    std::chrono::microseconds>(
                                                              end_time
                                                              - start_time)
                                                              .count();
    fprintf(out_fp, "# TIME\t%s (ms):\t%zu\n",msg,this->to_ms(elapsed_micro));
    start_time = end_time;
    return elapsed_micro;
  }
  size_t show_fp(FILE *out_fp, const std::string &msg)
  {
   return this->show_fp(out_fp, msg.c_str());
  }
  size_t show(const char *msg)
  {
    return show_fp(stdout,msg);
  }
  size_t show(const std::string &msg)
  {
    return show_fp(stdout,msg);
  }
  size_t locked_show(std::mutex *cout_mutex,const char *msg)
  {
    assert(cout_mutex != nullptr);
    const std::scoped_lock<std::mutex> cout_lock(*cout_mutex);
    return this->show(msg);
  }
  size_t locked_show(std::mutex *cout_mutex,const std::string &msg)
  {
    assert(cout_mutex != nullptr);
    const std::scoped_lock<std::mutex> cout_lock(*cout_mutex);
    return this->show(msg);
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
