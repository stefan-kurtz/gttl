#include <vector>
#include <thread>
#include <memory>
#include <iostream>
#include "utilities/unused.hpp"
#include "threading/thread_pool.hpp"
#include "threading/thread_pool_var.hpp"

static size_t fibonacci(size_t n)
{
  return n<=1 ? 1 : (fibonacci(n-1) + fibonacci(n-2));
}

static void sum_fibonacci_var(size_t thread_id, GTTL_UNUSED size_t task_id,
                              size_t n, size_t *thread_sums)
{
  thread_sums[thread_id] += fibonacci(n);
}

class SumFibThreadData
{
  private:
  size_t num_threads;
  public:
  size_t n;
  std::vector<size_t> thread_sums;
  SumFibThreadData(size_t _num_threads,size_t _n)
    : num_threads(_num_threads)
    , n(_n)
    , thread_sums(num_threads,0)
  {}
  void output(void) const noexcept
  {
    for (size_t idx = 0; idx < num_threads; idx++)
    {
      std::cout << "thread\t" << idx << "\t" << thread_sums[idx] << std::endl;
    }
  }
  ~SumFibThreadData(void) { }
};

static void sum_fibonacci(size_t thread_id, GTTL_UNUSED size_t task_id,
                          SumFibThreadData *thread_data)
{
  thread_data->thread_sums[thread_id] += fibonacci(thread_data->n);
}

int main(int argc, char *argv[])
{
  long readlong1, readlong2;

  if (argc != 3 || sscanf(argv[1], "%ld", &readlong1) != 1 || readlong1 < 0 ||
      sscanf(argv[2], "%ld", &readlong2) != 1 || readlong2 < 0)
  {
    std::cerr << "Usage: " << argv[0] << " <num_threads> <n>"
              << std::endl;
    return EXIT_FAILURE;
  }
  const size_t num_threads = (size_t) readlong1;
  const size_t n = (size_t) readlong2;
  std::vector<size_t> thread_sums(num_threads,0);
  GttlThreadPoolVar(num_threads,10,sum_fibonacci_var,n,thread_sums.data());
  for (size_t idx = 0; idx < num_threads; idx++)
  {
    std::cout << "thread\t" << idx << "\t" << thread_sums[idx] << std::endl;
  }
  SumFibThreadData thread_data(num_threads,n);
  GttlThreadPool(num_threads,10,sum_fibonacci,&thread_data);
  thread_data.output();
  return EXIT_SUCCESS;
}
