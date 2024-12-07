#include <iostream>
#include "sequences/gttl_fastq_iterator.hpp"
#include "utilities/gttl_line_iterator.hpp"
#include "iterator_compare.hpp"

int main()
{
  size_t pseudo_hash = 0;
  for(unsigned long long i = 0; i < ITERATIONS; i++)
  {
    constexpr const size_t buf = (1 << 14);
    GttlLineIterator<buf> line_it(FILENAME);
    GttlFastQIterator<GttlLineIterator<buf>> fastq_it(line_it);

    for (auto& q : fastq_it)
    {
      pseudo_hash += std::ranges::count(q.header_get(), '@');
      pseudo_hash += std::ranges::count(q.sequence_get(), 'a');
      pseudo_hash += std::ranges::count(q.quality_get(), 'F');
    }
  }
  std::cout << pseudo_hash << std::endl;
}
