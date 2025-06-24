/*
  Copyright (c) 2024 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2024 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cassert>
#include <string_view>
#include <cstdio>
#include <vector>
#include "sequences/gttl_fasta_generator.hpp"
#include "utilities/gttl_mmap.hpp"
#include "sequences/format_sequence.hpp"
#include "seq_reader_options.hpp"

template<class Iterator>
static void process_iterator(Iterator &iterator,
                             bool statistics,
                             size_t line_width)
{
  size_t seqnum = 0;
  size_t total_length = 0;
  for (auto &&si : iterator)
  {
    const std::string_view &sequence = si->sequence_get();
    if (line_width > 0)
    {
      const std::string_view &header = si->header_get();
      std::cout << '>' << header << '\n';
      gttl_format_sequence(sequence,line_width);
    }
    seqnum++;
    total_length += sequence.size();
  }
  if (statistics)
  {
    std::cout << "# number of sequences\t" << seqnum << '\n';
    std::cout << "# total length\t" << total_length << '\n';
    std::cout << "# mean length\t" << total_length / seqnum << '\n';
  }
}

int main(int argc,char *argv[])
{
  SeqReaderOptions options{0,false};

  try
  {
    options.parse(argc, argv);
  }
  catch (std::invalid_argument &e) /* check_err.py */
  {
    std::cerr << argv[0] << ": " << e.what() << '\n';
    return EXIT_FAILURE;
  }
  if (options.help_option_is_set())
  {
    return EXIT_SUCCESS;
  }
  try
  {
    const bool statistics = options.statistics_option_is_set();
    const size_t line_width = options.line_width_get();
    const std::vector<std::string> &inputfiles = options.inputfiles_get();
    constexpr const size_t buf_size = size_t{1} << size_t{14};
    for (const auto & inputfile : inputfiles)
    {
      if (options.mapped_option_is_set())
      {
        const Gttlmmap<char> mapped_file(inputfile.c_str());
        GttlFastAGenerator<buf_size> gttl_si(std::string_view(
                                               mapped_file.ptr(),
                                               mapped_file.size()));
        process_iterator<GttlFastAGenerator<buf_size>>(gttl_si,
                                                       statistics,line_width);
      } else
      {
        GttlFastAGenerator<buf_size> gttl_si(inputfile.c_str());
        process_iterator<GttlFastAGenerator<buf_size>>(gttl_si,statistics,
                                                       line_width);

      }
    }
  }
  catch (const std::exception &err)
  {
    std::cerr << argv[0] << ": " << err.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
