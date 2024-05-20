#ifndef FASTQ_PARTS_HPP
#define FASTQ_PARTS_HPP
#include <string>
#include <iostream>
#include <algorithm>
#include <string_view>
#include <stdexcept>
#include "utilities/gttl_mmap.hpp"

/* The follwing function was contributed by Florian Jochens and
   slightly modified by Stefan Kurtz. */

template<typename MappedBaseType>
static inline const MappedBaseType *fastq_next_read_start(
  const MappedBaseType *guess,
  const MappedBaseType *end_of_string)
{
  for (const MappedBaseType *ptr = guess; ptr < end_of_string; ptr++)
  {
    if (*ptr == '@')
    {
      const MappedBaseType *find_plus = ptr;
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

class FastQParts
{
  using MappedBaseType = char;
  Gttlmmap<MappedBaseType> mapped_file;
  const MappedBaseType *file_contents;
  std::vector<std::string_view> intervals;
  public:
  FastQParts(size_t num_parts,const std::string &inputfilename)
    : mapped_file(inputfilename.c_str())
    , file_contents(mapped_file.ptr())
    , intervals({})
  {
    if (inputfilename.size() >= 3 and
        inputfilename.substr(inputfilename.size() - 3) == std::string(".gz"))
    {
      throw std::string("cannot process .gz file with multiple threads");
    }
    assert(num_parts > 0 && mapped_file.size() > 0);
    if (num_parts > mapped_file.size())
    {
      intervals.push_back(std::string_view(file_contents, mapped_file.size()));
      return;
    }
    const size_t part_size = mapped_file.size()/num_parts;
    const MappedBaseType *end_of_mapped_string
      = file_contents + mapped_file.size();
    const MappedBaseType *guess = file_contents + part_size;
    size_t previous_start = 0;
    while(true)
    {
      const MappedBaseType *read_start
        = fastq_next_read_start<MappedBaseType>(guess,
                                                end_of_mapped_string);
      if (read_start == nullptr)
      {
        intervals.push_back(std::string_view(file_contents + previous_start,
                                             mapped_file.size() -
                                               previous_start));
        break;
      }
      const size_t this_start
        = static_cast<size_t>(read_start - file_contents);
      intervals.push_back(std::string_view(file_contents + previous_start,
                                           this_start - previous_start));
      previous_start = this_start;
      guess = read_start + part_size;
    }
  }
  void show(void) const noexcept
  {
    std::cout << "# fields: interval_start, interval_size" << std::endl;
    for (auto &&sw : intervals)
    {
      std::cout << static_cast<size_t>(sw.data() - file_contents)
                << "\t" << sw.size() << std::endl;
    }
  }
  using ConstIterator = std::vector<std::string_view>::const_iterator;
  ConstIterator begin(void) const
  {
    return intervals.cbegin();
  }

  ConstIterator end(void) const
  {
    return intervals.cend();
  }
  const std::string_view operator [](size_t idx) const
  {
    return intervals[idx];
  }
  size_t size(void) const
  {
    return intervals.size();
  }
  size_t total_size(void) const
  {
    return mapped_file.size();
  }
};
#endif
