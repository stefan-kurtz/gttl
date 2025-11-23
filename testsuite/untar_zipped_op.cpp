#include <stdexcept>
#include <string>
#include <iostream>
#include <vector>
#include "utilities/cxxopts.hpp"
#include "untar_zipped_op.hpp"

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << '\n';
}

UnzippedTarOptions::UnzippedTarOptions(void) = default;

void UnzippedTarOptions::parse(int argc, char **argv)
{
  cxxopts::Options options(argv[0],"unzip tar file and process content of "
                                   "files in memory; files must have suffix "
                                   ".tar.gz or .tar.bz2");
  options.set_width(80);
  options.custom_help(std::string("[options] inputfile0 [inputfile1...]"));
  options.set_tab_expansion();
  options.add_options()
    ("m,store_in_memory",
     "store the unzipped files as strings in memory",
     cxxopts::value<bool>(store_option)->default_value("false"))
    ("n,no_rapidgzip",
     "do not use rapidgzip, even if available",
     cxxopts::value<bool>(no_rapidgzip_option)->default_value("false"))
    ("h,help", "Print usage information");
  try
  {
    auto result = options.parse(argc, argv);
    if (result.count("help") > 0)
    {
      help_option = true;
      usage(options);
    }
    const std::vector<std::string>& unmatched_args = result.unmatched();
    if (unmatched_args.size() < 1)
    {
      throw std::invalid_argument("missing positional inputfiles "
                                  "(at least one filename is required)");
    }
    for (auto &&inputfile : unmatched_args)
    {
      inputfiles.push_back(inputfile);
    }
  }
  catch (const cxxopts::exceptions::exception &e)
  {
    usage(options);
    throw std::invalid_argument(e.what());
  }
}

const std::vector<std::string> &UnzippedTarOptions::inputfiles_get(void)
  const noexcept
{
  return inputfiles;
}

bool UnzippedTarOptions::store_option_is_set(void) const noexcept
{
  return store_option;
}

bool UnzippedTarOptions::help_option_is_set(void) const noexcept
{
  return help_option;
}

bool UnzippedTarOptions::no_rapidgzip_option_is_set(void) const noexcept
{
  return no_rapidgzip_option;
}
