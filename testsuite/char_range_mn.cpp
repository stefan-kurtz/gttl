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
#include <cassert>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <string>
#include <algorithm>
#include <iostream>
#include "utilities/cxxopts.hpp"
#include "utilities/mathsupport.hpp"
#include "sequences/gttl_seq_iterator.hpp"
#include "sequences/guess_if_protein_seq.hpp"
#include "sequences/char_range.hpp"

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << std::endl;
}

class CharRangeOptions
{
 private:
  std::vector<std::string> input_files{};
  bool invert_option = false, reverse_option = false, help_option = false;

 public:
  CharRangeOptions() {};

  void parse(int argc, char **argv)
  {
    cxxopts::Options options(argv[0],"compute ranges of characters of the same "
                                     "kind (like nucleotides or wildcard)");
    options.set_width(80);
    options.custom_help(std::string("[options] filename1 [filename1 ..]"));
    options.set_tab_expansion();
    options.add_options()
       ("i,invert", "show ranges of chararacters not in the given set",
        cxxopts::value<bool>(invert_option)->default_value("false"))
       ("r,reverse", "enumerate ranges in reverse order",
        cxxopts::value<bool>(reverse_option)->default_value("false"))
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
        this->input_files.push_back(unmatched_args[idx]);
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
  bool invert_option_is_set(void) const noexcept
  {
    return invert_option;
  }
  bool reverse_option_is_set(void) const noexcept
  {
    return reverse_option;
  }
};

template<bool forward,bool invert>
static void display_char_ranges(const char *inputfilename)
{
  constexpr const int buf_size = 1 << 14;
  const bool is_protein = guess_if_protein_file(inputfilename);

  if (is_protein)
  {
    throw std::string(": can only handle DNA sequences");
    /* check_err.py checked */
  }
  GttlFpType in_fp = gttl_fp_type_open(inputfilename,"rb");
  if (in_fp == nullptr)
  {
    throw std::string(": cannot open file");
    /* check_err.py checked */
  }
  GttlSeqIterator<buf_size> gttl_si(in_fp);
  size_t non_wildcard_ranges_total_length = 0;
  try /* need this, as the catch needs to close the file pointer
         to prevent a memory leak */
  {
    size_t seqnum = 0;
    for (auto &&si : gttl_si)
    {
      auto sequence = std::get<1>(si);
      static constexpr const char nucleotides[] = "ACGTacgt";
      GttlCharRange<forward,invert,nucleotides> ranger(sequence.data(),
                                                    sequence.size());
      for (auto &&it : ranger)
      {
        std::cout << seqnum << "\t" << std::get<0>(it)
                  << "\t" << std::get<1>(it) << std::endl;
        non_wildcard_ranges_total_length += std::get<1>(it);
      }
      seqnum++;
    }
  }
  catch (std::string &msg)
  {
    gttl_fp_type_close(in_fp);
    throw msg;
  }
  gttl_fp_type_close(in_fp);
  std::cout << "# non_wildcard_ranges_total_length\t" <<
                non_wildcard_ranges_total_length << std::endl;
}

int main(int argc,char *argv[])
{
  CharRangeOptions options;

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
  const char *progname = argv[0];
  bool haserr = false;
  if (options.invert_option_is_set())
  {
    for (int idx = 2; idx < argc; idx++)
    {
      try
      {
        display_char_ranges<true,true>(argv[idx]);
      }
      catch (std::string &msg)
      {
        std::cerr << progname << ": file \"" << argv[idx] << "\""
                  << msg << std::endl;
        haserr = true;
        break;
      }
    }
  } else
  {
    for (int idx = 1; idx < argc; idx++)
    {
      try
      {
        display_char_ranges<true,false>(argv[idx]);
      }
      catch (std::string &msg)
      {
        std::cerr << progname << ": file \"" << argv[idx] << "\""
                  << msg << std::endl;
        haserr = true;
        break;
      }
    }
  }
  return haserr ? EXIT_FAILURE : EXIT_SUCCESS;
}
