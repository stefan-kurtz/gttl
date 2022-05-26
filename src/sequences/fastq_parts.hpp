#ifndef GTTL_FASTQ_PARTS_HPP
#define GTTL_FASTQ_PARTS_HPP
#include <string>
#include <iostream>
#include <algorithm>
#include "utilities/gttl_mmap.hpp"

/* The follwing function was contributed by Florian Jochens and
   slightly modified by Stefan Kurtz. */

static inline const char *fastq_next_read_start(const char *guess,
                                                const char *end_of_string)
{
  for (const char *ptr = guess; ptr < end_of_string; ptr++)
  {
    if (*ptr == '@')
    {
      const char *find_plus = ptr;
      /* now check that the current line is the header line. This is verified
         by reading the next two occurrences of \n and checking that the
         next char is +, that is the next but one line is the third line of
         a fastq entry */
      /* find first EOL */
      for (/* Nothing */; find_plus < end_of_string && *find_plus != '\n';
           find_plus++)
           /* Nothing */ ;
      /* find next EOL */
      find_plus++;
      for (/* Nothing */; find_plus < end_of_string && *find_plus != '\n';
           find_plus++)
           /* Nothing */ ;
      if (find_plus + 1 < end_of_string && find_plus[1] == '+')
      {
        return ptr;
      }
    }
  }
  return nullptr;
}

struct FastQParts
{
  Gttlmmap<char> mapped_file;
  const char *file_contents;
  std::vector<std::pair<size_t,size_t>> intervals;
  FastQParts(size_t num_parts,const std::string &inputfilename)
    : mapped_file(Gttlmmap<char>(inputfilename.c_str()))
    , file_contents(mapped_file.ptr())
    , intervals({})
  {
    assert(num_parts > 0 && mapped_file.size() > 0);
    const size_t part_size = mapped_file.size()/num_parts;
    const char *end_of_mapped_string = file_contents + mapped_file.size();
    const char *guess = file_contents + part_size;
    size_t previous_start = 0;
    while(true)
    {
      const char *read_start = fastq_next_read_start(guess,
                                                     end_of_mapped_string);
      if (read_start == nullptr)
      {
        intervals.push_back({previous_start,
                             mapped_file.size() - previous_start});
        break;
      }
      const size_t this_start
        = static_cast<size_t>(read_start - file_contents);
      intervals.push_back({previous_start,this_start - previous_start});
      previous_start = this_start;
      guess = read_start + part_size;
    }
  }
  void show(void) const noexcept
  {
    std::cout << "# fields: interval_start, interval_end" << std::endl;
    for (auto &&itv : intervals)
    {
      std::cout << std::get<0>(itv) << "\t" << std::get<1>(itv) << std::endl;
    }
  }
};
#endif
