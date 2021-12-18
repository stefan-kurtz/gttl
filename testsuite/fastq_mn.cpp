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
#include "utilities/cxxopts.hpp"
#include "sequences/gttl_fastq_iterator.hpp"

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << std::endl;
}

class FastQReaderOptions
{
 private:
  std::vector<std::string> inputfiles{};
  bool help_option = false, statistics_option = false,
       echo_option = false, fasta_output_option = false;
  size_t split_size = 0;

 public:
  FastQReaderOptions(void) {};

  void parse(int argc, char **argv)
  {
    cxxopts::Options options(argv[0],"process fastq files, optionally output "
                                     "statistics, echo in the input or show "
                                     "sequences in fasta file");
    options.set_width(80);
    options.custom_help(std::string("[options] filename0 [filename1]"));
    options.set_tab_expansion();
    options.add_options()
       ("s,statistics", "output statistics",
        cxxopts::value<bool>(statistics_option)->default_value("false"))
       ("e,echo", "available only for single file",
        cxxopts::value<bool>(echo_option)->default_value("false"))
       ("split_size", "specify number of sequences for each split",
        cxxopts::value<size_t>(split_size)->default_value("0"))
       ("f,fasta_output", "output sequences in fasta format; when two files "
                          "are specified, the sequences are zipped, i.e. "
                          "sequences at even indexes (when counting from 0) "
                          "are from the first file and sequences at odd "
                          "indexes are from the second file",
        cxxopts::value<bool>(fasta_output_option)->default_value("false"))
       ("h,help", "print usage");
    try
    {
      auto result = options.parse(argc, argv);
      if (result.count("help") > 0)
      {
        help_option = true;
        usage(options);
      }
      const std::vector<std::string>& unmatched_args = result.unmatched();
      for (size_t idx = 0; idx < unmatched_args.size(); idx++)
      {
        inputfiles.push_back(unmatched_args[idx]);
      }
      if (inputfiles.size() < 1)
      {
        throw cxxopts::OptionException("not enough input files");
      }
      if (inputfiles.size() > 2)
      {
        throw cxxopts::OptionException("superfluous input files");
      }
      if (split_size != 0 && inputfiles.size() != 1)
      {
        throw cxxopts::OptionException("option --splitsize_t is only available "
                                       "for a single file");
      }
    }
    catch (const cxxopts::OptionException &e)
    {
      usage(options);
      throw std::invalid_argument(e.what());
    }
  }
  bool help_option_is_set(void) const noexcept
  {
    return help_option;
  }
  bool statistics_option_is_set(void) const noexcept
  {
    return statistics_option;
  }
  bool echo_option_is_set(void) const noexcept
  {
    return echo_option;
  }
  bool fasta_output_option_is_set(void) const noexcept
  {
    return fasta_output_option;
  }
  size_t split_size_get(void) const noexcept
  {
    return split_size;
  }
  const std::vector<std::string> &inputfiles_get(void) const noexcept
  {
    return inputfiles;
  }
};

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

static void process_single_file(bool statistics,
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

static void process_paired_files(bool statistics,
                                 bool fasta_output,
                                 const std::string &filename0,
                                 const std::string &filename1)
{
  constexpr const int buf_size = 1 << 14;
  GttlFastQIterator<buf_size> fastq_it0(filename0),
                              fastq_it1(filename1);

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
        process_single_file(statistics,echo,fasta_output,inputfiles[0]);
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
