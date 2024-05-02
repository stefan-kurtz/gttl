#ifndef RANDOM_SAMPLE_HPP
#define RANDOM_SAMPLE_HPP
#include <cstddef>
#include <cassert>
#include <algorithm>
#include <vector>
#include "utilities/uniform_random_int.hpp"
#include "utilities/multibitvector.hpp"

template<typename BaseType>
static inline std::vector<BaseType> gttl_random_sample(size_t total_number,
                                                       size_t sample_size,
                                                       unsigned int seed)
{
  assert(total_number > 0);
  UniformRandomInteger<BaseType> uri(0,total_number-1,seed);
  std::vector<BaseType> sample{};
  Multibitvector<false> occurs(total_number);
  while (sample.size() < sample_size)
  {
    const BaseType rseqnum = uri.get();
    if (not occurs[rseqnum])
    {
      sample.push_back(rseqnum);
      occurs.set(rseqnum);
    }
  }
  std::sort(sample.begin(),sample.end());
  return sample;
}
#endif
