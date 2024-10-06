#include <cstdint>
#include <stdexcept>
#include <string>
#include "utilities/cxxopts.hpp"
#include "ntcard_opt.hpp"

static void usage(const cxxopts::Options& options)
{
  std::cerr << options.help() << std::endl;
}

NtcardOptions::NtcardOptions(){};

void NtcardOptions::parse(int argc, char** argv)
{
  cxxopts::Options options(argv[0], "");
  options.set_width(80);
  options.custom_help(std::string("[options] inputfile"));
  options.set_tab_expansion();
  options.add_options()
     ("q,qgram_length", "specify qgram_length",
      cxxopts::value<size_t>(qgram_length)->default_value("32"))
     ("s", "specify s", cxxopts::value<size_t>(s)->default_value("7"))
     ("r", "specify r", cxxopts::value<size_t>(r)->default_value("27"))
     ("l,long", "Show histogram.",
      cxxopts::value<bool>(show_f_option)->default_value("false"))
     ("w,wildcard_as_A", "Handle a wildcard character just like character A",
      cxxopts::value<bool>(handle_wildcard_like_A)->default_value("false"))
     ("f,fast", "Only compute F0.",
      cxxopts::value<bool>(fast_option)->default_value("false"))
     ("b,binary", "Use binary NtTable.",
      cxxopts::value<bool>(binary_option)->default_value("false"))
     ("t,threads", "Number of threads to use.",
      cxxopts::value<size_t>(num_threads)->default_value("1"))
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
      help_option = false;
      const std::vector<std::string>& unmatched_args = result.unmatched();
      if (unmatched_args.size() < 1)
      {
        throw cxxopts::OptionException("missing input file");
      }
      if (unmatched_args.size() > 1)
      {
        throw cxxopts::OptionException("superfluous input file");
      }
      inputfile = unmatched_args[0];
    }
  }
  catch (const cxxopts::OptionException& e)
  {
    throw std::invalid_argument(e.what());
  }
}

bool NtcardOptions::help_option_is_set(void) const noexcept
{
  return help_option;
}

const std::string& NtcardOptions::inputfile_get(void) const noexcept
{
  return inputfile;
}

size_t NtcardOptions::qgram_length_get(void) const noexcept
{
  return qgram_length;
}

size_t NtcardOptions::s_get(void) const noexcept { return s; }

size_t NtcardOptions::r_get(void) const noexcept { return r; }

bool NtcardOptions::show_f_option_is_set(void) const noexcept
{
  return show_f_option;
}

bool NtcardOptions::fast_option_is_set(void) const noexcept
{
  return fast_option;
}

bool NtcardOptions::binary_option_is_set(void) const noexcept
{
  return binary_option;
}

bool NtcardOptions::handle_wildcard_like_A_is_set(void) const noexcept
{
  return handle_wildcard_like_A;
}

size_t NtcardOptions::num_threads_get(void) const noexcept
{
  return num_threads;
}
