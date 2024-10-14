#ifndef SPLIT_HPP
#define SPLIT_HPP
#include <cstring>
#include <string_view>
#include <stdexcept>
#include "utilities/has_suffix_or_prefix.hpp"
#include "utilities/gttl_mmap.hpp"

static inline size_t fastq_next_read_start(const char *file_contents,
                                           size_t total_size,
                                           size_t current)
{
  const char *end_of_string = file_contents + total_size;
  for (size_t idx = current; idx < total_size; idx++)
  {
    if (file_contents[idx] == '@')
    {
      const char *find_plus = file_contents + idx;
      /* now check that the current line is the header line. This is verified
         by reading the next two occurrences of \n and checking that the
         next char is +, that is the next but one line is the third line of
         a fastq entry */
      /* find first EOL */
      for (/* Nothing */; find_plus < end_of_string and *find_plus != '\n';
           find_plus++)
           /* Nothing */ ;
      /* find next EOL */
      find_plus++;
      for (/* Nothing */; find_plus < end_of_string && *find_plus != '\n';
           find_plus++)
           /* Nothing */ ;
      if (find_plus + 1 < end_of_string && find_plus[1] == '+')
      {
        return idx;
      }
    }
  }
  return total_size;
}

class SequencesSplit
{
  Gttlmmap<char> mapped_file;
  const char *file_contents;
  std::vector<std::string_view> intervals;

  public:
  SequencesSplit(size_t num_parts, const std::string &inputfilename,
                 bool fasta_format)
    : mapped_file(inputfilename.c_str())
    , file_contents(mapped_file.ptr())
    , intervals({})
  {
    if (gttl_has_suffix(inputfilename,std::string(".gz")))
    {
      throw std::string("cannot process .gz");
    }
    assert(num_parts > 0 && mapped_file.size() > 0);
    if (num_parts > mapped_file.size())
    {
      intervals.push_back(std::string_view(file_contents, mapped_file.size()));
      return;
    }
    const size_t part_size = mapped_file.size() / num_parts;
    size_t current_start = 0;
    for (size_t idx = 1; idx < num_parts and current_start < mapped_file.size();
         idx++)
    {
      size_t current = std::max(part_size * idx, current_start);
      if (fasta_format)
      {
        while (current < mapped_file.size() and file_contents[current] != '>')
        {
          current++;
        }
      } else
      {
        current = fastq_next_read_start(file_contents,mapped_file.size(),
                                        current);
      }
      assert(current_start < current);
      intervals.push_back(std::string_view(file_contents + current_start,
                                           current - current_start));
      current_start = current;
    }
    if (current_start < mapped_file.size())
    {
      intervals.push_back(std::string_view(file_contents + current_start,
                                           mapped_file.size() - current_start));
    }
  }
  double variance(void) const
  {
    double v = 0.0;
    const double mean = static_cast<double>(total_size())/intervals.size();
    for (auto &&sw : intervals)
    {
      const double diff = sw.size() - mean;
      v += (diff * diff);
    }
    return v;
  }
  void show(void) const noexcept
  {
    printf("# Fields: interval_start, interval_size\n");
    for (auto &&sw : intervals)
    {
      printf("%zu\t%zu\n",static_cast<size_t>(sw.data() - file_contents),
                          sw.size());
    }
    printf("# variance\t%.0e\n",variance());
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
