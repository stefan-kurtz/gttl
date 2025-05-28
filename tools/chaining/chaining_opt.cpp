#include "utilities/cxxopts.hpp"
#include "chaining_opt.hpp"
#include <iostream>
#include <string>
#include <vector>

ChainingOptions::ChainingOptions(void)
  : inputfiles({})
  , help_option(false)
  , local_option(false)
  , silent_option(false)
{}

/*
Usage: bin/gt chain2dim [options] -m matchfile
Chain pairwise matches.

-m       Specify file containing the matches
         mandatory option
         default: undefined
-global  perform global chaining
         - optional parameter gc switches
           on gap costs (according to L1-model)
         - optional parameter ov means
           that overlaps between matches are allowed
         - optional parameter all means
           that all optimal chains are processed
-local   perform local chaining
         compute local chains (according to L1-model).
         - If no parameter is given, compute local chains with
           maximums score.
         - If parameter is given, this must be a positive number
           optionally followed by the character b or p.
         - If only the number, say k, is given, this is the minimum
           score of the chains output.
         - If a number is followed by character b, then output all
           chains with the largest k scores.
         - If a number is followed by character p, then output all
           chains with scores at most k percent away
           from the best score.
-wf      specify weight factor > 0.0 to obtain score of a fragment
         requires one of the options
         -local const
         -global gc
         -global ov
         default: 1.00
-maxgap  specify maximal width of gap in chain
         default: 0
-silent  do not output the chains but only report their lengths and scores
         default: no
-v       be verbose
         default: no
-help    display help and exit
-version display version information and exit
*/

void ChainingOptions::parse(int argc, char **argv)
{
  cxxopts::Options options(argv[0], "Chain pairwise matches");
  options.set_width(80);
  options.custom_help(std::string("[options] matchfile"));
  options.set_tab_expansion();
  options.add_options()
    ("l,local", "perform local chaining instead of global",
     cxxopts::value<bool>(local_option)->default_value("false"))
    ("s,silent", "do not output the chains but only report their lengths and "
                 "scores",
     cxxopts::value<bool>(silent_option)->default_value("false"))
    ("h,help", "print usage");
  try
  {
    auto result = options.parse(argc, argv);
    if (result.count("help") > 0)
    {
      help_option = true;
    }

    const std::vector<std::string>& unmatched_args = result.unmatched();
    for (const auto & unmatched_arg : unmatched_args)
    {
      inputfiles.push_back(unmatched_arg);
    }
    if (inputfiles.size() < 1)
    {
      throw cxxopts::exceptions::exception("not enough inputfiles");
    }
    if (inputfiles.size() > 2)
    {
      throw cxxopts::exceptions::exception("superfluous inputfiles");
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
bool ChainingOptions::help_option_is_set(void) const noexcept
{
  return help_option;
}
bool ChainingOptions::local_option_is_set(void) const noexcept
{
  return local_option;
}
bool ChainingOptions::silent_option_is_set(void) const noexcept
{
  return silent_option;
}
const std::string &ChainingOptions::inputfile_get(void)
  const noexcept
{
  return inputfiles[0];
}
