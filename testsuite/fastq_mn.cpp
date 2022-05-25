/*
  Copyright (c) 2021 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2021 Center for Bioinformatics, University of Hamburg

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
#include <iostream>
#include <fstream>
#include <cstdbool>
#include <typeinfo>
#include <cassert>
#include <algorithm>
#include <cstdio>
#include "utilities/str_format.hpp"
#include "utilities/basename.hpp"
#include "sequences/gttl_fastq_iterator.hpp"
#include "fastq_opt.hpp"

static void fastq_split_writer(size_t split_size,
                               const std::string &inputfilename)
{
  constexpr const int buf_size = 1 << 14;
  GttlFastQIterator<buf_size> fastq_it(inputfilename);
  int file_number = 0;
  auto it = fastq_it.begin();
  bool exhausted = false;
  GttlBasename inputfile_base(inputfilename.c_str());
  while (!exhausted)
  {
    if (it == fastq_it.end())
    {
      break;
    }
    StrFormat outfilename("%s_%02d",inputfile_base.str(),file_number++);
    std::ofstream out_stream;
    out_stream.open(outfilename.str());
    for (size_t idx = 0; idx < split_size; idx++)
    {
      if (it == fastq_it.end())
      {
        exhausted = true;
        break;
      }
      const std::string_view &sequence = (*it).sequence_get();
      const std::string_view &header = (*it).header_get();
      const std::string_view &quality = (*it).quality_get();
      out_stream << header << std::endl
                 << sequence << std::endl
                 << "+" << std::endl
                 << quality << std::endl;
      ++it;
    }
    out_stream.close();
  }
}

static void process_single_file_streamed(bool statistics,
                                         bool echo,
                                         bool fasta_output,
                                         const std::string &inputfilename)
{
  constexpr const int buf_size = 1 << 14;
  GttlFastQIterator<buf_size> fastq_it(inputfilename);
  size_t seqnum = 0, total_length = 0;

  for (auto &&fastq_entry : fastq_it)
  {
    const std::string_view &sequence = fastq_entry.sequence_get();
    if (echo)
    {
      const std::string_view &header = fastq_entry.header_get();
      const std::string_view &quality = fastq_entry.quality_get();
      std::cout << header << std::endl
                << sequence << std::endl
                << "+" << std::endl
                << quality << std::endl;
    } else
    {
      if (fasta_output)
      {
        const std::string_view &header = fastq_entry.header_get();
        std::cout << ">" << header.substr(1) << std::endl
                  << sequence << std::endl;
      }
    }
    total_length += sequence.size();
    seqnum++;
  }
  if (statistics)
  {
    std::cout << "# number of sequences\t" << seqnum << std::endl;
    std::cout << "# total length\t" << total_length << std::endl;
    std::cout << "# mean length\t" << total_length/seqnum << std::endl;
  }
}

#ifdef NEWCODE
static void process_single_file_mapped(bool statistics,
                                       bool echo,
                                       bool fasta_output,
                                       const std::string &inputfilename)
{
  Gttlmmap<char> mapped_file(filename);
  assert(mapped_file.size() > 0);
  const char *file_contents = mapped_file.ptr();
  GttlFastQIterator<buf_size> fastq_it(file_contents,mapped_file.size());
  size_t seqnum = 0, total_length = 0;

  for (auto &&fastq_entry : fastq_it)
  {
    const ModifiableStringView &sequence = fastq_entry.sequence_get();
    if (echo)
    {
      const std::string_view &header = fastq_entry.header_get();
      const std::string_view &quality = fastq_entry.quality_get();
      std::cout << header << std::endl
                << sequence << std::endl
                << "+" << std::endl
                << quality << std::endl;
    } else
    {
      if (fasta_output)
      {
        const std::string_view &header = fastq_entry.header_get();
        std::cout << ">" << header.substr(1) << std::endl
                  << sequence << std::endl;
      }
    }
    total_length += sequence.size();
    seqnum++;
  }
  if (statistics)
  {
    std::cout << "# number of sequences\t" << seqnum << std::endl;
    std::cout << "# total length\t" << total_length << std::endl;
    std::cout << "# mean length\t" << total_length/seqnum << std::endl;
  }
}
#endif

static void process_paired_files(bool statistics,
                                 bool fasta_output,
                                 const std::string &filename0,
                                 const std::string &filename1)
{
  constexpr const int buf_size = 1 << 14;
  GttlFastQIterator<buf_size> fastq_it0(filename0),fastq_it1(filename1);

  size_t seqnum = 0, total_length[2] = {0};
  auto it0 = fastq_it0.begin();
  auto it1 = fastq_it1.begin();
  while (it0 != fastq_it0.end() && it1 != fastq_it1.end())
  {
    const std::string_view &sequence0 = (*it0).sequence_get();
    const std::string_view &sequence1 = (*it1).sequence_get();
    if (fasta_output)
    {
      const std::string_view &header0 = (*it0).header_get();
      const std::string_view &header1 = (*it1).header_get();
      std::cout << ">" << header0.substr(1) << std::endl;
      std::cout << sequence0 << std::endl;
      std::cout << ">" << header1.substr(1) << std::endl;
      std::cout << sequence1 << std::endl;
    }
    total_length[0] += sequence0.size();
    total_length[1] += sequence1.size();
    seqnum++;
    ++it0;
    ++it1;
  }
  if (statistics)
  {
    std::cout << "# number of sequences\t" << seqnum << std::endl;
    std::cout << "# total length0\t" << total_length[0] << std::endl;
    std::cout << "# total length1\t" << total_length[1] << std::endl;
  }
}

int main(int argc,char *argv[])
{
  FastQReaderOptions options{};

  try
  {
    options.parse(argc, argv);
  }
  catch (std::invalid_argument &e) /* check_err.py */
  {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  if (options.help_option_is_set())
  {
    return EXIT_SUCCESS;
  }
  try
  {
    const bool statistics = options.statistics_option_is_set();
    const bool echo = options.echo_option_is_set();
    const bool fasta_output = options.fasta_output_option_is_set();
    const size_t split_size = options.split_size_get();
    const std::vector<std::string> &inputfiles = options.inputfiles_get();
    if (inputfiles.size() == 1)
    {
      if (split_size > 0)
      {
        fastq_split_writer(split_size,inputfiles[0]);
      } else
      {
        process_single_file_streamed(statistics,echo,fasta_output,
                                     inputfiles[0]);
      }
    } else
    {
      process_paired_files(statistics,fasta_output,inputfiles[0],inputfiles[1]);
    }
  }
  catch (std::string &msg)
  {
    std::cerr << argv[0] << ": " << msg << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
