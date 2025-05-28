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
  , window_size(1)
  , number_of_threads(1)
  , hash_bits(-1)
  , canonical_option(false)
  , at_constant_distance_option(false)
  , sort_by_hash_value_option(false)
  , help_option(false)
{}

void MinimizerOptions::parse(int argc, char **argv)
{
  const std::string default_window_size = "1",
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
    ("w,window_size", "number of qgrams in window",
     cxxopts::value<size_t>(window_size)
                 ->default_value(default_window_size))
    ("t,number_of_threads", "specify number of threads",
     cxxopts::value<size_t>(number_of_threads)->default_value("1"))
    ("b,hash_bits", "specify number of bits used for hashing, if undefined "
                   "(i.e. -1), then it is set to 2 * kmer_length",
     cxxopts::value<int>(hash_bits)->default_value("-1"))
    ("m,show_mode", "specify if or how to show the minimizer: 0 means no show, "
                    "1 means to use .show() of HashedQgram-Class and "
                    "2 means to use the iterator of the HashedQgram Class",
     cxxopts::value<int>(show_mode)->default_value("0"))
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
        throw cxxopts::exceptions::exception("not enough input files");
      }
      if (hash_bits == -1)
      {
        hash_bits = 2 * qgram_length;
      } else
      {
        if (hash_bits < static_cast<int>(2 * qgram_length))
        {
          StrFormat msg("hash_bits = %d < %zu = 2 * kmer_length is "
                        "not possible",hash_bits,2 * qgram_length);
          throw cxxopts::exceptions::exception(msg.str());
        }
      }
      if (show_mode != 0 and show_mode != 1 and show_mode != 2)
      {
        throw cxxopts::exceptions::exception(std::string("option -m,--show_mode must "
                                                   "be used with argument "
                                                   "0, 1  or 2"));
      }
    }
  }
  catch (const cxxopts::exceptions::exception &e)
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
  return window_size;
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

int MinimizerOptions::show_mode_get(void) const noexcept
{
  return show_mode;
}

bool MinimizerOptions::help_option_is_set(void) const noexcept
{
  return help_option;
}
