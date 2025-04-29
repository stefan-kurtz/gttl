#include <iostream>
#include <vector>
#include "utilities/cxxopts.hpp"
#include "utilities/str_format.hpp"
#include "sequences/multiseq_factory.hpp"

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << std::endl;
}

#define CHECK_PAIRWISE_EXCLUDE(I_IDX,J_IDX)\
        if (std::get<1>(values[I_IDX]) > 0 and std::get<1>(values[J_IDX]) > 0)\
        {\
          StrFormat msg("option %s and %s exclude each other", \
                        std::get<0>(values[I_IDX]),std::get<0>(values[J_IDX]));\
          throw std::invalid_argument(msg.str());\
        }

class MultiseqFactoryOptions
{
  private:
  std::vector<std::string> inputfiles;
  size_t num_parts,
         len_parts,
         num_sequences,
         sequence_output_width;
  bool statistics_option,
       help_option;

  public:
  MultiseqFactoryOptions(void)
    : num_parts(0)
    , len_parts(0)
    , num_sequences(0)
    , sequence_output_width(0)
    , statistics_option(false)
    , help_option(false)
  {}
  void parse(int argc, char **argv)
  {
    cxxopts::Options options(argv[0],"create MultiSeqFactory from inputfiles; "
                                     "if two files are given as argument, then "
                                     "it is supposed that they are paired end"
                                     " fastq files");
    options.set_width(80);
    options.custom_help(std::string("[options] inputfile0 [inputfile1]"));
    options.set_tab_expansion();
    options.add_options()
      ("p,num_parts",
       "Number of parts to split the file into "
       "(mutually exclusive with -l and -n)",
       cxxopts::value<size_t>(num_parts)->default_value("0"))
      ("l,len_parts",
       "Length of parts to split the file into "
       "(mutually exclusive with -p and -n)",
       cxxopts::value<size_t>(len_parts)->default_value("0"))
      ("n,num_sequences",
       "Number of sequences per file (mutually exclusive with -l and -p)",
       cxxopts::value<size_t>(num_sequences)->default_value("0"))
      ("s,statistics", "output statistics of the sizes of the different parts",
        cxxopts::value<bool>(statistics_option)->default_value("false"))
      ("w,width", "output sequences in lines of width specified by the "
                  "argument of this option",
        cxxopts::value<size_t>(sequence_output_width)->default_value("0"))
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
                                    "(one or two are required)");
      } else
      {
        if (unmatched_args.size() > 2)
        {
          throw std::invalid_argument("superfluous positional inputfiles"
                                      "(one or two are required)");
        }
      }
      for (auto &&inputfile : unmatched_args)
      {
        inputfiles.push_back(inputfile);
      }
      const std::vector<std::pair<const char *,size_t>> values{
        std::make_pair("-p",num_parts),
        std::make_pair("-l",len_parts),
        std::make_pair("-n",num_sequences)};
      CHECK_PAIRWISE_EXCLUDE(0,1);
      CHECK_PAIRWISE_EXCLUDE(0,2);
      CHECK_PAIRWISE_EXCLUDE(1,2);
    }
    catch (const cxxopts::OptionException &e)
    {
      usage(options);
      throw std::invalid_argument(e.what());
    }
    if (inputfiles.size() == 2 and num_sequences > 0 and
        num_sequences % size_t(2) != 0)
    {
      throw std::invalid_argument("if option -n is used for two files, the "
                                  "arguments must be an even number");
    }
  }
  size_t num_parts_get(void) const noexcept
  {
    return num_parts;
  }
  size_t len_parts_get(void) const noexcept
  {
    return len_parts;
  }
  size_t num_sequences_get(void) const noexcept
  {
    return num_sequences;
  }
  size_t sequence_output_width_get(void) const noexcept
  {
    return sequence_output_width;
  }
  const std::vector<std::string> &inputfiles_get(void) const noexcept
  {
    return inputfiles;
  }
  bool help_option_is_set(void) const noexcept
  {
    return help_option;
  }
  bool statistics_option_is_set(void) const noexcept
  {
    return statistics_option;
  }
};

static void test_multiseq_factory(size_t num_parts,
                                  size_t len_parts,
                                  size_t num_sequences,
                                  size_t sequence_output_width,
                                  bool statistics_option,
                                  const std::vector<std::string> &inputfiles)
{
  const uint8_t padding_char = UINT8_MAX;
  const bool short_header = true;
  GttlMultiseqFactory *multiseq_factory = nullptr;
  if (inputfiles.size() == 1)
  {
    multiseq_factory
      = new GttlMultiseqFactory(inputfiles[0],
                                num_parts,
                                len_parts,
                                num_sequences,
                                padding_char,
                                short_header);
  } else
  {
    assert (inputfiles.size() == 2);
    multiseq_factory
      = new GttlMultiseqFactory(inputfiles[0],
                                inputfiles[1],
                                num_parts,
                                len_parts,
                                num_sequences,
                                padding_char,
                                short_header);
  }
  std::cout << "# number of parts\t" << multiseq_factory->size() << std::endl;
  if (statistics_option)
  {
    multiseq_factory->statistics();
  }
  if (sequence_output_width > 0)
  {
    multiseq_factory->sequence_output(sequence_output_width);
  }
  delete multiseq_factory;
}

int main(int argc, char *argv[])
{
  MultiseqFactoryOptions options;
  try
  {
    options.parse(argc,argv);
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
  try
  {
    test_multiseq_factory(options.num_parts_get(),
                          options.len_parts_get(),
                          options.num_sequences_get(),
                          options.sequence_output_width_get(),
                          options.statistics_option_is_set(),
                          options.inputfiles_get());
  }
  catch (std::string &msg)
  {
    std::cerr << argv[0] << msg << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
