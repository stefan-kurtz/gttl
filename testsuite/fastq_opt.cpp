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
#include <cassert>
#include <algorithm>
#include "utilities/cxxopts.hpp"
#include "fastq_opt.hpp"

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << std::endl;
}

FastQReaderOptions::FastQReaderOptions(void) {};

void FastQReaderOptions::parse(int argc, char **argv)
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
     ("m,mapped", "use fastq reader which works on mapped file",
      cxxopts::value<bool>(mapped_option)->default_value("false"))
     ("split_size", "specify number of sequences for each split",
      cxxopts::value<size_t>(split_size)->default_value("0"))
     ("t,threads", "specify number of threads for parallel reading",
      cxxopts::value<size_t>(num_threads)->default_value("0"))
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
      throw cxxopts::OptionException("option --splitsize is only available "
                                     "for a single file");
    }
    if (num_threads != 0)
    {
      if (split_size != 0)
      {
        throw cxxopts::OptionException("option -t/--threads and --splitsize "
                                       "are not compatible");
      }
      if (echo_option)
      {
        throw cxxopts::OptionException("option -t/--threads and --echo "
                                       "are not compatible");
      }
      if (fasta_output_option)
      {
        throw cxxopts::OptionException("option -t/--threads and "
                                       "-f/--fasta_output are not compatible");
      }
    }
  }
  catch (const cxxopts::OptionException &e)
  {
    usage(options);
    throw std::invalid_argument(e.what());
  }
}

bool FastQReaderOptions::help_option_is_set(void) const noexcept
{
  return help_option;
}

bool FastQReaderOptions::statistics_option_is_set(void) const noexcept
{
  return statistics_option;
}

bool FastQReaderOptions::echo_option_is_set(void) const noexcept
{
  return echo_option;
}

bool FastQReaderOptions::fasta_output_option_is_set(void) const noexcept
{
  return fasta_output_option;
}

bool FastQReaderOptions::mapped_option_is_set(void) const noexcept
{
  return mapped_option;
}

size_t FastQReaderOptions::split_size_get(void) const noexcept
{
  return split_size;
}

size_t FastQReaderOptions::num_threads_get(void) const noexcept
{
  return num_threads;
}

const std::vector<std::string> &FastQReaderOptions::inputfiles_get(void)
      const noexcept
{
  return inputfiles;
}
