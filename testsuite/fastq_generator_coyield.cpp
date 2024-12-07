#include "sequences/gttl_fastq_generator.hpp"
#include <algorithm>
#include <cassert>
#include <iostream>
#include "iterator_compare.hpp"


int main()
{
  size_t pseudo_hash = 0;
  for(unsigned long long i = 0; i < ITERATIONS; i++)
  {
    for (const FastQEntry && q : gttl_read_fastq(FILENAME))
    {
      pseudo_hash += std::ranges::count(q.header_get(), '@');
      pseudo_hash += std::ranges::count(q.sequence_get(), 'a');
      pseudo_hash += std::ranges::count(q.quality_get(), 'F');
    }
  }
  std::cout << pseudo_hash << std::endl;
}
