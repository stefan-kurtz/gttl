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
#include <cstdbool>
#include <cassert>
#include <algorithm>
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

 public:
  FastQReaderOptions(void) {};

  void parse(int argc, char **argv)
  {
    cxxopts::Options options(argv[0],"process fastaq files, optionally output "
                                     "statistics, echo in the input or show "
                                     "sequences in fasta file");
    options.set_width(80);
    options.custom_help(std::string("[options] filename0 [filename1]"));
    options.set_tab_expansion();
    options.add_options()
       ("s,statistics", "output statistics",
        cxxopts::value<bool>(statistics_option)->default_value("false"))
       ("e,echo", "available only for file",
        cxxopts::value<bool>(echo_option)->default_value("false"))
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
  const std::vector<std::string> &inputfiles_get(void) const noexcept
  {
    return inputfiles;
  }
};

static void process_single_file(bool statistics,
                                bool echo,
                                bool fasta_output,
                                const std::string &filename)
{
  constexpr const int buf_size = 1 << 14;
  GttlFastQIterator<buf_size> fastq_it(filename);
  size_t seqnum = 0, total_length = 0;

  for (auto &&fastq_entry : fastq_it)
  {
    const std::string &sequence = std::get<1>(fastq_entry);
    if (echo)
    {
      const std::string &quality = std::get<3>(fastq_entry);
      std::cout << fastq_header(fastq_entry) << std::endl;
      std::cout << sequence << std::endl;
      std::cout << "+" << std::endl;
      std::cout << quality << std::endl;
    }
    if (fasta_output)
    {
      std::cout << ">" << seqnum << std::endl;
      std::cout << sequence << std::endl;
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
    const std::string &sequence0 = fastq_sequence(*it0);
    const std::string &sequence1 = fastq_sequence(*it1);
    if (fasta_output)
    {
      std::cout << ">" << seqnum << ".0" << std::endl;
      std::cout << sequence0 << std::endl;
      std::cout << ">" << seqnum << ".1" << std::endl;
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
    bool statistics = options.statistics_option_is_set();
    bool echo = options.echo_option_is_set();
    bool fasta_output = options.fasta_output_option_is_set();
    const std::vector<std::string> &inputfiles = options.inputfiles_get();
    if (inputfiles.size() == 1)
    {
      process_single_file(statistics,echo,fasta_output,inputfiles[0]);
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
