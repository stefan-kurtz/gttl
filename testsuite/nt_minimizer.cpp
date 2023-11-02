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

#include "utilities/cxxopts.hpp"
#include "sequences/complement_plain.hpp"
#include "sequences/qgrams_hash_nthash.hpp"
#include "sequences/gttl_multiseq.hpp"
#include "sequences/qgrams_hash_nthash.hpp"
#include "sequences/hashed_qgrams.hpp"

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << std::endl;
}

class NtMinimizerOptions
{
 private:
  std::vector<std::string> inputfiles;
  size_t qgram_length,
         minimum_mem_length,
         number_of_threads;
  int hashbits;
  bool canonical_option,
       at_constant_distance_option,
       sort_by_hash_value_option,
       help_option;

 public:
  NtMinimizerOptions(void)
    : inputfiles({})
    , qgram_length(0)
    , minimum_mem_length(0)
    , number_of_threads(1)
    , hashbits(-1)
    , canonical_option(false)
    , at_constant_distance_option(false)
    , sort_by_hash_value_option(false)
    , help_option(false)
  {};

  void parse(int argc, char **argv)
  {
    const std::string default_minimum_mem_length = "30",
                      default_qgram_length = "18";

    cxxopts::Options options(argv[0],"generate minimizers of DNA sequences "
                                     "in given input files");
    options.set_width(80);
    options.custom_help(std::string("[options] filename0 [filename1]"));
    options.set_tab_expansion();
    options.add_options()
      ("k,kmer_length", "specify k-mer length",
       cxxopts::value<size_t>(qgram_length)
                      ->default_value(default_qgram_length))
      ("l,minimum_mem_length", "minimum length of maximal exact match (MEM)",
       cxxopts::value<size_t>(minimum_mem_length)
                   ->default_value(default_minimum_mem_length))
      ("t,number_of_threads", "specify number of threads",
       cxxopts::value<size_t>(number_of_threads)->default_value("1"))
      ("b,hashbits", "specify number of bits used for hashing, if undefined "
                     "(i.e. -1), then it is set to 2 * kmer_length",
       cxxopts::value<int>(hashbits)->default_value("-1"))
      ("c,canonical", "use canonical minimizers",
       cxxopts::value<bool>(canonical_option)->default_value("false"))
      ("d,constant_distance", "compute hashed_qgrams at constant distance",
       cxxopts::value<bool>(at_constant_distance_option)
                ->default_value("false"))
      ("s,sort_by_hash_value", "sort the hashed qgrams in ascending order of "
                               "their hash value",
       cxxopts::value<bool>(sort_by_hash_value_option)
                ->default_value("false"))
     ("h,help", "print usage");
    try
    {
      auto result = options.parse(argc, argv);
      if (result.count("help") > 0)
      {
        help_option = true;
        usage(options);
      } else
      {
        const std::vector<std::string>& unmatched_args = result.unmatched();
        for (size_t idx = 0; idx < unmatched_args.size(); idx++)
        {
          inputfiles.push_back(unmatched_args[idx]);
        }
        if (inputfiles.size() < 1)
        {
          throw cxxopts::OptionException("not enough input files");
        }
        if (hashbits == -1)
        {
          hashbits = 2 * qgram_length;
        } else
        {
          if (hashbits < static_cast<int>(2 * qgram_length))
          {
            StrFormat msg("hashbits = %d < %lu = 2 * kmer_length is "
                          "not possible",hashbits,2 * qgram_length);
            throw cxxopts::OptionException(msg.str());
          }
        }
        if (minimum_mem_length < qgram_length)
        {
          StrFormat msg("minimum_mem_length = %lu, but it must be larger "
                        "than kmer_length = %lu",
                        minimum_mem_length,qgram_length);
          throw cxxopts::OptionException(msg.str());
        }
      }
    }
    catch (const cxxopts::OptionException &e)
    {
      usage(options);
      throw std::invalid_argument(e.what());
    }
  }
  const std::vector<std::string> &inputfiles_get(void) const noexcept
  {
    return inputfiles;
  }
  size_t qgram_length_get(void) const noexcept
  {
    return qgram_length;
  }
  size_t window_size_get(void) const noexcept
  {
    assert(minimum_mem_length >= qgram_length_get());
    return minimum_mem_length - qgram_length_get() + 1;
  }
  size_t number_of_threads_get(void) const noexcept
  {
    return number_of_threads;
  }
  int hashbits_get(void) const noexcept
  {
    return hashbits;
  }
  bool canonical_option_is_set(void) const noexcept
  {
    return canonical_option;
  }
  bool at_constant_distance_option_is_set(void) const noexcept
  {
    return at_constant_distance_option;
  }
  bool sort_by_hash_value_option_is_set(void) const noexcept
  {
    return sort_by_hash_value_option;;
  }
  bool help_option_is_set(void) const noexcept
  {
    return help_option;
  }
};

void run_nt_minimizer(const NtMinimizerOptions &options)
{
  RunTimeClass rt_create_multiseq{};
  GttlMultiseq *multiseq = nullptr;
  try
  {
    constexpr const bool store_sequences = true;
    constexpr const uint8_t padding_char = UINT8_MAX;
    if (options.canonical_option_is_set())
    {
      multiseq = multiseq_with_reverse_complement<store_sequences,
                                                  complement_plain>
                                                 (options.inputfiles_get(),
                                                  padding_char);
    } else
    {
      multiseq = new GttlMultiseq(options.inputfiles_get(),
                                  store_sequences,
                                  padding_char);
    }
  }
  catch (std::string &msg) /* check_err.py */
  {
    delete multiseq;
    throw msg;
  }
  rt_create_multiseq.show("reading input files and creating multiseq");
  int hashbits;
  int var_sizeof_unit_hashed_qgram;
  assert(options.hashbits_get() != -1);
  for (auto &log : multiseq->statistics())
  {
    std::cout << "# " << log << std::endl;
  }
  if (multiseq->sequences_bits_get() + options.hashbits_get() <= 64)
  {
    hashbits = 64 - multiseq->sequences_bits_get();
    var_sizeof_unit_hashed_qgram = 8;
  } else
  {
    if (multiseq->sequences_bits_get() + options.hashbits_get() <= 72)
    {
      hashbits = 72 - multiseq->sequences_bits_get();
      var_sizeof_unit_hashed_qgram = 9;
    } else
    {
      StrFormat msg(": file \"%s\", cannot handle "
                    "hashbits + sequence_bits = %d + %d > 72",
                    options.inputfiles_get()[0].c_str(),
                    options.hashbits_get(),
                    multiseq->sequences_bits_get());
      delete multiseq;
      throw msg;
    }
  }
  assert(var_sizeof_unit_hashed_qgram == 8 or
         var_sizeof_unit_hashed_qgram == 9);
  constexpr_for<8,9+1,1>([&](auto sizeof_unit_hashed_qgram)
  {
    if (sizeof_unit_hashed_qgram == var_sizeof_unit_hashed_qgram)
    {
      std::vector<std::string> log_vector{};
      if (options.canonical_option_is_set())
      {
        using HashedQgrams = HashedQgramsGeneric<sizeof_unit_hashed_qgram,
                                                 QgramNtHashIterator4>;

        HashedQgrams hqg (*multiseq,
                          options.number_of_threads_get(),
                          options.qgram_length_get(),
                          options.window_size_get(),
                          hashbits,
                          options.sort_by_hash_value_option_is_set(),
                          options.at_constant_distance_option_is_set(),
                          &log_vector);
        RunTimeClass rt_output_hashed_qgrams{};
        hqg.show();
        log_vector.push_back(rt_output_hashed_qgrams
                             .to_string("output of hashed kmers"));
      } else
      {
        using HashedQgrams = HashedQgramsGeneric<sizeof_unit_hashed_qgram,
                                                 QgramNtHashFwdIterator4>;
        HashedQgrams hqg (*multiseq,
                          options.number_of_threads_get(),
                          options.qgram_length_get(),
                          options.window_size_get(),
                          hashbits,
                          options.sort_by_hash_value_option_is_set(),
                          options.at_constant_distance_option_is_set(),
                          &log_vector);
        RunTimeClass rt_output_hashed_qgrams{};
        hqg.show();
        log_vector.push_back(rt_output_hashed_qgrams
                             .to_string("output of hashed kmers"));
      }
      for (auto &msg : log_vector)
      {
        std::cout << "# " << msg << std::endl;
      }
    }
  });
  delete multiseq;
}

int main(int argc, char *argv[])
{
  RunTimeClass rt_total{};
  NtMinimizerOptions options{};

  try
  {
    options.parse(argc, argv);
  }
  catch (std::invalid_argument &e) /* check_err.py */
  {
    std::cerr << argv[0] << ": file " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  if (options.help_option_is_set())
  {
    return EXIT_SUCCESS;
  }
  RunTimeClass rt_create_multiseq{};
  try
  {
    run_nt_minimizer(options);
  }
  catch (const std::string &msg)
  {
    for (auto &&inputfile : options.inputfiles_get())
    {
      std::cerr << argv[0] << ": file \"" << inputfile << "\""
                << msg << std::endl;
    }
    return EXIT_FAILURE;
  }
  rt_total.show("overall");
  return EXIT_SUCCESS;
}
