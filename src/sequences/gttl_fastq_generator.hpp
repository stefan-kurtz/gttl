#ifndef GTTL_FASTQ_GENERATOR_HPP
#define GTTL_FASTQ_GENERATOR_HPP

#include <generator>
#include "utilities/gttl_line_generator_coyield.hpp"

struct FastQEntry
{
  // Constructor using move-semantics
  FastQEntry(std::string&& h, std::string&& s, std::string&& q)
      : header(std::move(h))
      , sequence(std::move(s))
      , quality(std::move(q))
  {}

  // We return only a view, since that is almost always more efficient.
  std::string_view header_get() const noexcept
  {
    return std::string_view(header);
  }
  std::string_view sequence_get() const noexcept
  {
    return std::string_view(sequence);
  }
  std::string_view quality_get() const noexcept
  {
    return std::string_view(quality);
  }

 private:
  // We store the data as std::string, thus holding ownership in the struct.
  std::string header;
  std::string sequence;
  std::string quality;
};

template <const size_t buf_size = (1 << 14)>
std::generator<FastQEntry> gttl_read_fastq(GttlFpType fp)
{
  size_t state = 0;
  // We need to use strings here, since a view into line would be overwritten
  // on each iteration of the loop.
  // However, this is fine because we can std::move() the strings themselves.
  std::string header, sequence, quality;

  for (const auto&& line : gttl_read_lines<buf_size>(fp))
  {
    switch (state)
    {
      case 0:
        header = line;
        break;
      case 1:
        sequence = line;
        break;
      case 2:
        break;
      case 3:
        quality = line;
        // Construct the entry using move-semantics
        co_yield FastQEntry(std::move(header),
                            std::move(sequence),
                            std::move(quality));
        state = 0;
        continue;
    }
    state++;
  }
}

std::generator<FastQEntry> gttl_read_fastq(const std::string file_path)
{
  return gttl_read_fastq(gttl_fp_type_open(file_path.c_str(), "r"));
}

#endif  // GTTL_FASTQ_GENERATOR_HPP
