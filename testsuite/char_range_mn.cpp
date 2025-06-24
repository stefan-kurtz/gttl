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
#include <cstdint>
#include <exception>
#include <stdexcept>
#include <string>
#include <iostream>
#include <vector>
#include "sequences/alphabet.hpp"
#include "sequences/gttl_fasta_generator.hpp"
#include "sequences/gttl_multiseq.hpp"
#include "utilities/cxxopts.hpp"
#include "sequences/guess_if_protein_seq.hpp"
#include "sequences/char_range.hpp"
#include "sequences/char_finder.hpp"
#include "sequences/literate_multiseq.hpp"
#include "utilities/gttl_file_open.hpp"

#ifdef _WIN32
  #include "utilities/windows_getopt.hpp"
#else
  // We disable the include-cleaner check here and when using optind, since
  // optind is defind in <bits/getopt_core.h>, which should *NOT* be
  // included directly.
  #include <unistd.h> // NOLINT(misc-include-cleaner)
#endif

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << std::endl;
}

class CharRangeOptions
{
 private:
  std::vector<std::string> inputfiles{};
  bool help_option = false, multiseq_option = false, singlechar_option = false,
       invert_option = false, reverse_option = false;

 public:
  CharRangeOptions() {};

  void parse(int argc, char **argv)
  {
    cxxopts::Options options(argv[0],"compute ranges of characters of the same "
                                     "kind (like nucleotides or wildcard)");
    options.set_width(80);
    options.custom_help(std::string("[options] filename1 [filename2 ...]"));
    options.set_tab_expansion();
    options.add_options()
       ("s,singlechar", "use single character finder for N",
        cxxopts::value<bool>(singlechar_option)->default_value("false"))
       ("i,invert", "show ranges of chararacters not in the given set",
        cxxopts::value<bool>(invert_option)->default_value("false"))
       ("r,reverse", "enumerate ranges in reverse order",
        cxxopts::value<bool>(reverse_option)->default_value("false"))
       ("m,multiseq", "create Multiseq from each file and process the entire "
                      "concatenated sequence as a single unit",
        cxxopts::value<bool>(multiseq_option)->default_value("false"))
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
        this->inputfiles.push_back(unmatched_args[idx]);
      }
    }
    catch (const cxxopts::exceptions::exception &e)
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
  bool singlechar_option_is_set(void) const noexcept
  {
    return singlechar_option;
  }
  bool multiseq_option_is_set(void) const noexcept
  {
    return multiseq_option;
  }
  const std::vector<std::string> &inputfiles_get(void) const noexcept
  {
    return this->inputfiles;
  }
};

template<class CharFinder,const CharFinder &char_finder,bool forward,
         bool invert>
static void display_char_ranges(const char *inputfilename)
{
  const bool is_protein = guess_if_protein_file(inputfilename);

  if (is_protein)
  {
    throw std::runtime_error(": can only handle DNA sequences");
    /* check_err.py checked */
  }
  GttlFpType in_fp = gttl_fp_type_open(inputfilename,"rb");
  if (in_fp == nullptr)
  {
    throw std::runtime_error(": cannot open file");
    /* check_err.py checked */
  }
  GttlFastAGenerator fasta_gen(in_fp);
  size_t ranges_total_length = 0;
  using ThisCharRange = GttlCharRange<CharFinder,char_finder,forward,invert>;
  size_t seqnum = 0;
  for (const auto *si : fasta_gen)
  {
    auto sequence = si->sequence_get();
    ThisCharRange ranger(sequence.data(),sequence.size());
    for (auto const &&range : ranger)
    {
      std::cout << seqnum << "\t" << std::get<0>(range)
                << "\t" << std::get<1>(range) << std::endl;
      ranges_total_length += std::get<1>(range);
    }
    seqnum++;
  }
  std::cout << "# ranges_total_length\t" << ranges_total_length << std::endl;
}

template<class CharFinder,const CharFinder &char_finder>
static bool display_char_ranges_cases(const char *progname,
                                      const CharRangeOptions &options)
{
  bool haserr = false;
  for (auto && inputfile : options.inputfiles_get())
  {
    std::cout << inputfile << std::endl;
    try
    {
      if (options.invert_option_is_set())
      {
        if (options.reverse_option_is_set())
        {
          display_char_ranges<CharFinder,char_finder,false,true>
                             (inputfile.c_str());
        } else
        {
          display_char_ranges<CharFinder,char_finder,true,true>
                             (inputfile.c_str());
        }
      } else
      {
        if (options.reverse_option_is_set())
        {
          display_char_ranges<CharFinder,char_finder,false,false>
                             (inputfile.c_str());
        } else
        {
          display_char_ranges<CharFinder,char_finder,true,false>
                             (inputfile.c_str());
        }
      }
    }
    catch (const std::exception &err)
    {
      std::cerr << progname << ": file \"" << inputfile << "\""
                << err.what() << std::endl;
      haserr = true;
      break;
    }
  }
  return haserr;
}

template<class CharFinder,const CharFinder &char_finder,
         bool forward,bool invert>
static void display_char_ranges_multiseq(const GttlMultiseq *encoded_multiseq)
{
  const char *sequence_ptr = encoded_multiseq->sequence_ptr_get();
  assert(encoded_multiseq->sequences_number_get() > 0);
  const size_t sequence_length_concatenation
    = encoded_multiseq->sequences_total_length_get() +
      encoded_multiseq->sequences_number_get() - 1;
  using ThisCharRange
    = GttlCharRange<CharFinder,char_finder,forward,invert>;
  ThisCharRange ranger(sequence_ptr,sequence_length_concatenation);
  size_t ranges_total_length = 0;
  for (auto &&range : ranger)
  {
    std::cout << std::get<0>(range)
              << "\t" << std::get<1>(range) << std::endl;
    ranges_total_length += std::get<1>(range);
  }
  std::cout << "# ranges_total_length\t" << ranges_total_length
            << std::endl;
}

template<class CharFinder,const CharFinder &char_finder>
static void display_char_ranges_multiseq_cases(const CharRangeOptions &options,
                                               GttlMultiseq *multiseq)
{
  if (options.invert_option_is_set())
  {
    if (options.reverse_option_is_set())
    {
      display_char_ranges_multiseq<CharFinder,char_finder,false,true>(multiseq);
    } else
    {
      display_char_ranges_multiseq<CharFinder,char_finder,true,true>(multiseq);
    }
  } else
  {
    if (options.reverse_option_is_set())
    {
      display_char_ranges_multiseq<CharFinder,char_finder,false,false>
                                  (multiseq);
    } else
    {
      display_char_ranges_multiseq<CharFinder,char_finder,true,false>
                                  (multiseq);
    }
  }
}

/* These objects are only used in certain blocks, but some
   compilers, like g++ (SUSE Linux) 7.5.0 require that they
   are in the global scope */

static constexpr const
  char_finder::EncodedNucleotideFinder encoded_nucleotide_finder;
static constexpr const char_finder::NFinder single_N_finder{};
static constexpr const char_finder::NucleotideFinder nucleotide_finder{};

int main(int argc,char *argv[])
{
  CharRangeOptions options;

  try
  {
    options.parse(argc, argv);
  }
  catch (const std::invalid_argument &err) /* check_err.py */
  {
    std::cerr << argv[0] << ": " << err.what() << std::endl;
    return EXIT_FAILURE;
  }
  if (options.help_option_is_set())
  {
    return EXIT_SUCCESS;
  }
  bool haserr = false;
  if (options.multiseq_option_is_set())
  {
    for (auto && inputfile : options.inputfiles_get())
    {
      std::cout << inputfile << std::endl;
      GttlMultiseq *multiseq = nullptr;
      try
      {
        constexpr const bool store_header = true;
        constexpr const bool store_sequence = true;
        multiseq = new GttlMultiseq(inputfile.c_str(),
                                    store_header,store_sequence,UINT8_MAX);
      }
      catch (const std::exception &err)
      {
        std::cerr << argv[0]
                  << ": file \""
                  << argv[optind] //NOLINT(misc-include-cleaner)
                  << "\""
                  << err.what()
                  << std::endl;
        haserr = true;
      }
      if (!haserr)
      {
        LiterateMultiseq<alphabet::nucleotides_upper_lower,4>
          lit_multiseq(multiseq);
        lit_multiseq.perform_sequence_encoding();

        display_char_ranges_multiseq_cases<char_finder::EncodedNucleotideFinder,
                                           encoded_nucleotide_finder>
                                          (options,multiseq);
      }
      delete multiseq;
    }
  } else
  {
    if (options.singlechar_option_is_set())
    {
      haserr = display_char_ranges_cases<char_finder::NFinder,single_N_finder>
                                        (argv[0],options);
    } else
    {
      haserr = display_char_ranges_cases<char_finder::NucleotideFinder,
                                         nucleotide_finder>
                                        (argv[0],options);
    }
  }
  return haserr ? EXIT_FAILURE : EXIT_SUCCESS;
}
