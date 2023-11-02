#include <cstddef>
#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>
#include "utilities/str_format.hpp"
#include "utilities/cxxopts.hpp"
#include "minimizer_opt.hpp"

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << std::endl;
}

MinimizerOptions::MinimizerOptions(void)
  : inputfiles({})
  , qgram_length(0)
  , minimum_mem_length(0)
  , number_of_threads(1)
  , hash_bits(-1)
  , canonical_option(false)
  , at_constant_distance_option(false)
  , sort_by_hash_value_option(false)
  , help_option(false)
{}

void MinimizerOptions::parse(int argc, char **argv)
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
    ("b,hash_bits", "specify number of bits used for hashing, if undefined "
                   "(i.e. -1), then it is set to 2 * kmer_length",
     cxxopts::value<int>(hash_bits)->default_value("-1"))
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
      if (hash_bits == -1)
      {
        hash_bits = 2 * qgram_length;
      } else
      {
        if (hash_bits < static_cast<int>(2 * qgram_length))
        {
          StrFormat msg("hash_bits = %d < %lu = 2 * kmer_length is "
                        "not possible",hash_bits,2 * qgram_length);
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

const std::vector<std::string> &MinimizerOptions::inputfiles_get(void)
  const noexcept
{
  return inputfiles;
}

size_t MinimizerOptions::qgram_length_get(void) const noexcept
{
  return qgram_length;
}

size_t MinimizerOptions::window_size_get(void) const noexcept
{
  assert(minimum_mem_length >= qgram_length_get());
  return minimum_mem_length - qgram_length_get() + 1;
}

size_t MinimizerOptions::number_of_threads_get(void) const noexcept
{
  return number_of_threads;
}

int MinimizerOptions::hash_bits_get(void) const noexcept
{
  return hash_bits;
}

bool MinimizerOptions::canonical_option_is_set(void) const noexcept
{
  return canonical_option;
}

bool MinimizerOptions::at_constant_distance_option_is_set(void) const noexcept
{
  return at_constant_distance_option;
}

bool MinimizerOptions::sort_by_hash_value_option_is_set(void) const noexcept
{
  return sort_by_hash_value_option;;
}

bool MinimizerOptions::help_option_is_set(void) const noexcept
{
  return help_option;
}
