#include <cstdio>
#include <algorithm>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cassert>
#include "utilities/fs_prio_store.hpp"
#include "utilities/runtime_class.hpp"

template<bool debug,bool heap_based>
static void run_trial(size_t num_elements,size_t max_value,size_t capacity)
{
  std::vector<int> all_elements{};
  FSPrioStore<int,heap_based> fspq(capacity);
  for(size_t i = 0; i < num_elements; i++)
  {
    if (debug)
    {
      fspq.show();
    }
    int v = rand() % (max_value+1);
    if (debug)
    {
      printf("new=%d\n",v);
    }
    all_elements.push_back(v);
    fspq.add(v);
  }
  std::sort(all_elements.begin(),all_elements.end());
  fspq.sort();
  size_t check_idx = all_elements.size();
  for (auto &v : fspq)
  {
    assert (check_idx > 0);
    check_idx--;
    if (all_elements[check_idx] != v)
    {
      std::cerr << "all_elements[" << check_idx << " = "
                << all_elements[check_idx] << " != " << v << std::endl;
      exit(EXIT_FAILURE);
    }
    if (debug)
    {
      std::cout << v << std::endl;
    }
  }
}

template<bool debug,bool heap_based>
static void run_trial_iter(size_t trials, size_t num_elements,
                           size_t max_value,size_t capacity)
{
  for (size_t t = 0; t < trials; t++)
  {
    run_trial<debug,heap_based>(num_elements,max_value, capacity);
  }
}

int main(int argc,char *argv[])
{
  int trials, num_elements, max_value, capacity;
  if (argc != 5 or
      std::sscanf(argv[1],"%d",&trials) != 1 or trials < 1 or
      std::sscanf(argv[2],"%d",&num_elements) != 1 or num_elements < 1 or
      std::sscanf(argv[3],"%d",&max_value) != 1 or max_value < 1 or
      std::sscanf(argv[4],"%d",&capacity) != 1 or capacity < 1)
  {
    fprintf(stderr,"Usage: %s <trials> <num_elements> <max_value> <capacity>\n",
            argv[0]);
    return EXIT_FAILURE;
  }
  constexpr const bool debug = false;
  RunTimeClass rt0{};
  run_trial_iter<debug,true>(static_cast<size_t>(trials),
                             static_cast<size_t>(num_elements),
                             static_cast<size_t>(max_value),
                             static_cast<size_t>(capacity));
  rt0.show("run_trials heap_based=true");
  RunTimeClass rt1{};
  run_trial_iter<debug,false>(static_cast<size_t>(trials),
                              static_cast<size_t>(num_elements),
                              static_cast<size_t>(max_value),
                              static_cast<size_t>(capacity));
  rt1.show("run_trials heap_based=false");
  return EXIT_SUCCESS;
}
