#include <vector>
#include <thread>
#include <memory>
#include <iostream>
#include "utilities/unused.hpp"
#include "threading/thread_pool_var.hpp"

void test(size_t thread_id, GTTL_UNUSED size_t task_id,
          size_t num_elements, size_t *thread_sums)
{
  size_t this_sum = 0;
  for (size_t i = 0; i < num_elements; i++)
  {
    this_sum += i;
  }
  thread_sums[thread_id] = this_sum;
}

int main(int argc, char *argv[])
{
  long readlong1, readlong2;

  if (argc != 3 || sscanf(argv[1], "%ld", &readlong1) != 1 || readlong1 < 0 ||
      sscanf(argv[2], "%ld", &readlong2) != 1 || readlong2 < 0)
  {
    std::cerr << "Usage: " << argv[0] << " <num_threads> <elements>"
              << std::endl;
    return EXIT_FAILURE;
  }
  size_t num_threads = (size_t) readlong1;
  size_t num_elements = (size_t) readlong2;
  size_t *thread_sums = new size_t [num_threads];
  GttlThreadPool(num_threads,1000,test,num_elements,thread_sums);
  for (size_t idx = 0; idx < num_threads; idx++)
  {
    std::cout << "thread\t" << idx << "\t" << thread_sums[idx] << std::endl;
  }
  delete[] thread_sums;
  return EXIT_SUCCESS;
}
