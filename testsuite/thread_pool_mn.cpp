#include <vector>
#include <thread>
#include <memory>
#include <cstring>
#include <iostream>
#include "utilities/unused.hpp"
#include "threading/thread_pool_var.hpp"

static size_t fibonacci(size_t n)
{
  return n<=1 ? 1 : (fibonacci(n-1) + fibonacci(n-2));
}

static void sum_fibonacci(size_t thread_id, GTTL_UNUSED size_t task_id,
                          size_t n, size_t *thread_sums)
{
  thread_sums[thread_id] += fibonacci(n);
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
  size_t *thread_sums = new size_t [num_threads];
  memset(thread_sums,0,num_threads * sizeof *thread_sums);
  GttlThreadPool(num_threads,10,sum_fibonacci,n,thread_sums);
  for (size_t idx = 0; idx < num_threads; idx++)
  {
    std::cout << "thread\t" << idx << "\t" << thread_sums[idx] << std::endl;
  }
  delete[] thread_sums;
  return EXIT_SUCCESS;
}
