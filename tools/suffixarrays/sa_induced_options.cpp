/*
  Copyright (c) 2013-2025 Stefan Kurtz <stefan.kurtz@uni-hamburg.de>
  Copyright (c) 2013-2025 Center for Bioinformatics, University of Hamburg

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

#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <stdexcept>
#include <filesystem>
#include <vector>
#include <string>

#include "utilities/cxxopts.hpp"
#include "utilities/string_values_join.hpp"
#include "utilities/option_choices.hpp"
#include "utilities/split_string.hpp"
#include "sa_induced_options.hpp"

class DisplayOptions
{
  const std::vector<std::string> allowed_options{"abs_suftab",
                                                 "rel_suftab",
                                                 "lcptab",
                                                 "tistab"};
  std::vector<std::string> args_set;
  public:
  DisplayOptions(void) = default;

  [[nodiscard]] std::string allowed_help(void) const
  {
    return string_values_join(", ", allowed_options.begin(),
                                    allowed_options.end());
  }

  void set(const std::string &args)
  {
    if (args == std::string(""))
    {
      return;
    }
    const std::vector<std::string> split_args = gttl_split_string(args,',');
    for (auto &&arg : split_args)
    {
      if (std::ranges::find(allowed_options, arg)
            != allowed_options.end())
      {
        args_set.push_back(arg);
      } else
      {
        throw std::invalid_argument(std::string("illegal argument \"" + arg +
                                                "\" to option -s/--show, " +
                                                "possible values: " +
                                                allowed_help()));
      }
    }
  }
#define DISPLAY_FUNCTION(KEY)\
bool KEY##_output(void) const\
{\
  return std::ranges::find(args_set, #KEY) != args_set.end();\
}
  DISPLAY_FUNCTION(abs_suftab)
  DISPLAY_FUNCTION(rel_suftab)
  DISPLAY_FUNCTION(lcptab)
  DISPLAY_FUNCTION(tistab)
};

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << '\n';
}

void SainOptions::parse(int argc, char **argv)
{
  const std::map<std::string,std::string> lcptab_method_choices_map
      {{"","no lcptab construction and output"},
       {"kasai13n","use algorithm of kasai with 13n bytes space peak"},
       {"kasai9n","use semi external implementation of algorithm of kasai "
                  "with 9n bytes space peak"},
       {"plcp5n","use semi external implementation of algorithm involving "
                 "PLCP table with 5n bytes space peak"}};

  const OptionChoices lcptab_method_choices(lcptab_method_choices_map);
  DisplayOptions display_options{};
  cxxopts::Options options(
      argv[0],
      "construct suffix array using the linear time induced suffix "
      "sorting algorithm\n  of Yuta Mori, also implemented in the "
      "sais-lite software\n\n");
  options.set_width(80);
  options.custom_help(std::string("[options] filename1 [filename2 ...]"));
  options.set_tab_expansion();
  options.add_options()
     ("v,verbose", "run in verbose mode (default: no) ",
      cxxopts::value<bool>(verbose_opt)->default_value("false"))

     ("plain_input_format",
      "process input file as plain without format "
      "considerations (default: no)",
      cxxopts::value<bool>(plain_input_format_opt)
               ->default_value("false"))

     ("c,check_suftab", "check the generated suftab",
      cxxopts::value<bool>(check_suftab_opt)->default_value("false"))

     ("i,intset_sizeof",
      "set the size of the basetype of the representation of the inset "
      "to 0, 1 , 2 or 4; 0 means that the basetype which leads to minimum size "
      "is automatically determined; the other values stand for the number of "
      "bytes used for the base type (1 means uint8_t, 2 means uint16_t, "
      "4 means unt32_t",
      cxxopts::value<int>(intset_sizeof)->default_value("-1"))

     ("s,show",
      "comma separated string specifing what to show, possible values: " +
        display_options.allowed_help(),
      cxxopts::value<std::string>(show_options_spec)->default_value(""))

     ("a,absolute_suftab", "write the suftab with absolute positions to file",
      cxxopts::value<bool>(abs_suftab_out_opt)->default_value("false"))

     ("r,relative_suftab",
      "write the suftab with relative positions to file",
      cxxopts::value<bool>(rel_suftab_out_opt)->default_value("false"))

     ("l,lcptab",
      std::string("write the lcptab to file and specify the method to use; the "
                  "following methods are available: ")
                  + lcptab_method_choices.help_line(),
      cxxopts::value<std::string>(lcptab_method_string)
                                ->default_value(""))

     ("t,tistab", "write translated sequence to file",
      cxxopts::value<bool>(tistab_out_opt)->default_value("false"))

     ("reverse_complement", "compute index for all sequencees including "
      "their reverse complement",
      cxxopts::value<bool>(reverse_complement_option) ->default_value("false"))

     ("buffered", "use buffering when inserting Sstar suffixes",
      cxxopts::value<bool>(buffered_option)->default_value("false"))

     ("succinct", "create file with succinct representatin of lcp-table "
                  "instead of the byte representation with "
                  "saturated values stored in files with suffixes "
                  ".ll2 and .ll4",
      cxxopts::value<bool>(succinct_option)->default_value("false"))

     ("o,indexname",
      "set the name of the output files.\n"
      "option is mandatory, if more then one input file is given.",
      cxxopts::value<std::string>(indexname))

     ("h,help", "print usage");
  try
  {
    auto result = options.parse(argc, argv);
    if (result.contains("help"))
    {
      help_opt = true;
      usage(options);
    }
    const std::vector<std::string> &unmatched_args = result.unmatched();
    if (unmatched_args.empty())
    {
      throw std::invalid_argument("missing positional reference file argument");
    }
    for (const auto & unmatched_arg : unmatched_args)
    {
      inputfiles.push_back(unmatched_arg);
    }
    display_options.set(show_options_spec);
    abs_suftab_show_opt = display_options.abs_suftab_output();
    rel_suftab_show_opt = display_options.rel_suftab_output();
    tistab_show_opt = display_options.tistab_output();
    lcptab_show_opt = display_options.lcptab_output();
  }
  catch (const cxxopts::exceptions::exception &err)
  {
    usage(options);
    throw std::invalid_argument(err.what());
  }

  if (intset_sizeof != -1 && intset_sizeof != 0 &&
      intset_sizeof != 1 && intset_sizeof != 2 && intset_sizeof != 4)
  {
    throw std::invalid_argument("please select one of the values 0, 1, 2 or 4 "
                                "for option -i/--intset_sizeof.");
  }
  if (intset_sizeof != -1 && plain_input_format_opt)
  {
    std::cout << (" option -i/--intset_sizeof has no effect if option "
                  "--plain_input_format is set.")
              << '\n';
  }

  if (inputfiles.size() > 1 && indexname.empty())
  {
    throw std::invalid_argument("please select an indexname for the output"
                                " files via option -o/--indexname.");
  } else
  {
    if (inputfiles.size() == 1 && indexname.empty())
    {
      indexname = std::filesystem::path(inputfiles[0]).filename().string();
      assert(not indexname.empty());
    }
  }

  const int ret = lcptab_method_choices.choose(lcptab_method_string);
  if (ret == -1)
  {
    throw std::invalid_argument(std::string("illegal argument ") +
                                lcptab_method_string +
                                std::string("for option -l,--lcptab"));
  }
  lcptab_method = static_cast<LcptabMethod>(ret);
  if ((lcptab_method == Lcptab_kasai9n or lcptab_method == Lcptab_plcp5n)
      and not abs_suftab_out_opt_is_set())
  {
    throw std::invalid_argument("option -l/--lcptab kasai9n|plcp5n require "
                                "to use option -a/--absolute_suftab");
  }
  if (check_suftab_opt_is_set() and
      (lcptab_method == Lcptab_kasai9n or lcptab_method == Lcptab_plcp5n))
  {
    throw std::invalid_argument("option -l/--lcptab kasai9n|plcp5n and "
                                "option -c/--check_suftab are not compatible");
  }
  if ((not rel_suftab_out_opt) and rel_suftab_show_opt)
  {
    throw std::invalid_argument("display of relative_suftab requires to use "
                                "option -r/--relatative_suftab");
  }
  if (rel_suftab_out_opt and plain_input_format_opt)
  {
    throw std::invalid_argument("option -r/--relatative_suftab and "
                                "--plain_input_format are not compatible");
  }
  if (succinct_option and lcptab_method != Lcptab_plcp5n)
  {
    throw std::invalid_argument("option --succinct requires to use option "
                                "-l/--lcptab with argument plcp5n");
  }
}

bool SainOptions::help_opt_is_set(void) const noexcept
{
  return help_opt;
}

bool SainOptions::verbose_opt_is_set(void) const noexcept
{
  return verbose_opt;
}

bool SainOptions::plain_input_format_opt_is_set(void) const noexcept
{
  return plain_input_format_opt;
}

bool SainOptions::check_suftab_opt_is_set(void) const noexcept
{
  return check_suftab_opt;
}

bool SainOptions::abs_suftab_out_opt_is_set(void) const noexcept
{
  return abs_suftab_out_opt;
}

bool SainOptions::abs_suftab_show_opt_is_set(void) const noexcept
{
  return abs_suftab_show_opt;
}

bool SainOptions::rel_suftab_out_opt_is_set(void) const noexcept
{
  return rel_suftab_out_opt;
}

bool SainOptions::rel_suftab_show_opt_is_set(void) const noexcept
{
  return rel_suftab_show_opt;
}

bool SainOptions::lcptab_show_opt_is_set(void) const noexcept
{
  return lcptab_show_opt;
}

SainOptions::LcptabMethod SainOptions::lcptab_method_get(void) const noexcept
{
  return lcptab_method;
}

bool SainOptions::tistab_out_opt_is_set(void) const noexcept
{
  return tistab_out_opt;
}

bool SainOptions::tistab_show_opt_is_set(void) const noexcept
{
  return tistab_show_opt;
}

bool SainOptions::reverse_complement_option_is_set(void) const noexcept
{
  return reverse_complement_option;
}

bool SainOptions::buffered_option_is_set(void) const noexcept
{
  return buffered_option;
}

bool SainOptions::succinct_option_is_set(void) const noexcept
{
  return succinct_option;
}

std::string SainOptions::indexname_get(void) const noexcept
{
  return indexname;
}

const std::vector<std::string> &SainOptions::inputfiles_get(void) const noexcept
{
  return inputfiles;
}

int SainOptions::intset_sizeof_get(void) const noexcept
{
  return intset_sizeof;
}
