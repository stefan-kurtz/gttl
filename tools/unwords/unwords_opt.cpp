#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "utilities/cxxopts.hpp"
#include "unwords_opt.hpp"

UnwordsOptions::UnwordsOptions(void)
  : inputfiles({})
  , qgram_length_max(0)
  , help_option(false)
  , ignore_reverse_complement_option(false)
  , store_sequences_option(false)
{}

void UnwordsOptions::parse(int argc, char **argv)
{
  cxxopts::Options options(argv[0], "process fasta files and "
                                    "compute unwords");
  options.set_width(80);
  options.custom_help(std::string("[options] filename [filename...]"));
  options.set_tab_expansion();
  options.add_options()
    ("q,qgram_length_max", "specify maximal length of qgrams",
     cxxopts::value<size_t>(qgram_length_max)->default_value("0"))

    ("i,ignore_rc", "do not consider the reverse complement of the sequences",
     cxxopts::value<bool>(ignore_reverse_complement_option)
                         ->default_value("false"))

    ("s,store_sequences", "store sequences instead of reading them again "
                          "and again",
     cxxopts::value<bool>(store_sequences_option)->default_value("false"))

    ("h,help", "print usage");
  try
  {
    auto result = options.parse(argc, argv);
    if (result.contains("help"))
    {
      help_option = true;
    }

    const std::vector<std::string>& unmatched_args = result.unmatched();

    for (const auto & unmatched_arg : unmatched_args)
    {
      inputfiles.push_back(unmatched_arg);
    }
    if (inputfiles.empty())
    {
      throw cxxopts::exceptions::exception("not enough inputfiles");
    }
  }
  catch (const cxxopts::exceptions::exception &e)
  {
    std::cerr << options.help() << '\n';
    if (!help_option)
    {
      throw std::invalid_argument(e.what());
    }
  }
}

bool UnwordsOptions::help_option_is_set(void) const noexcept
{
  return help_option;
}

bool UnwordsOptions::ignore_reverse_complement_option_is_set(void)
  const noexcept
{
  return ignore_reverse_complement_option;
}

bool UnwordsOptions::store_sequences_option_is_set(void) const noexcept
{
  return store_sequences_option;
}

size_t UnwordsOptions::qgram_length_max_get(void) const noexcept
{
  return qgram_length_max;
}

const std::vector<std::string> &UnwordsOptions::inputfiles_get(void)
  const noexcept
{
  return inputfiles;
}
