#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include "utilities/cxxopts.hpp"
#include "utilities/runtime_class.hpp"
#include "sequences/literate_multiseq.hpp"
#include "utilities/random_sample.hpp"

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << std::endl;
}

class MultiseqOptions
{
 private:
  std::vector<std::string> inputfiles;
  bool help_option,
       protein_option,
       zipped_option,
       rankdist_option,
       short_header_option,
       statistics_option;
  size_t sample_size;
  unsigned int seed;
  int width_arg = -1;

 public:
  MultiseqOptions(void)
   : inputfiles({})
   , help_option(false)
   , protein_option(false)
   , zipped_option(false)
   , rankdist_option(false)
   , short_header_option(false)
   , statistics_option(false)
   , sample_size(0)
   , seed(0)
 {}

  void parse(int argc, char **argv)
  {
    cxxopts::Options options(argv[0],"tests for GttlMultiseq");
    options.set_width(80);
    options.custom_help(std::string("[options] filename1 [filename2 ...]"));
    options.set_tab_expansion();
    options.add_options()
       ("p,protein", "handle protein sequences",
        cxxopts::value<bool>(protein_option)->default_value("false"))
       ("sample", "extract random sample of the specified number of sequences",
        cxxopts::value<size_t>(sample_size)->default_value("0"))
       ("seed", "specify seed for generating random sample (default is 0, so "
                "that each call uses a differnt seed)",
        cxxopts::value<unsigned int>(seed)->default_value("0"))
       ("z,zipped", "expect two fastq  files with the same "
                    "number of sequences; show them "
                    "in zipped order, i.e. the "
                    "sequences at even indexes (when "
                    "counting from 0) are from the first "
                    "file and sequences at odd indexes are "
                    "from the second file",
        cxxopts::value<bool>(zipped_option)->default_value("false"))
       ("r,rankdist", "output distribution of ranks of "
                      "transformed sequences",
        cxxopts::value<bool>(rankdist_option)->default_value("false"))
       ("statistics", "output statistics about the sequences",
        cxxopts::value<bool>(statistics_option)->default_value("false"))
       ("s,short_header", "show header up to and excluding the first blank",
        cxxopts::value<bool>(short_header_option)->default_value("false"))
       ("w,width", "output headers and sequences; "
                   "width specifies the linewidth of the"
                   "sequence output; 0 means to output "
                   "a sequence in a single line",
        cxxopts::value<int>(width_arg)->default_value("-1"))
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
        if (unmatched_args.size() == 0)
        {
          throw std::invalid_argument("at least one inputput file is required");
        }
        for (size_t idx = 0; idx < unmatched_args.size(); idx++)
        {
          inputfiles.push_back(unmatched_args[idx]);
        }
      }
      if (zipped_option && inputfiles.size() != 2)
      {
        throw std::invalid_argument("option -z/--zipped requires exactly "
                                    "two files");
      }
    }
    catch (const cxxopts::OptionException &e)
    {
      usage(options);
      throw std::invalid_argument(e.what());
    }
  }
  bool help_option_is_set(void) const noexcept
  {
    return help_option;
  }
  bool protein_option_is_set(void) const noexcept
  {
    return protein_option;
  }
  bool zipped_option_is_set(void) const noexcept
  {
    return zipped_option;
  }
  bool rankdist_option_is_set(void) const noexcept
  {
    return rankdist_option;
  }
  bool short_header_option_is_set(void) const noexcept
  {
    return short_header_option;
  }
  bool statistics_option_is_set(void) const noexcept
  {
    return statistics_option;
  }
  int width_option_get(void) const noexcept
  {
    return width_arg;
  }
  const std::vector<std::string> &inputfiles_get(void) const noexcept
  {
    return inputfiles;
  }
  size_t sample_size_get(void) const noexcept
  {
    return sample_size;
  }
  unsigned int seed_get(void) const noexcept
  {
    return seed;
  }
};

int main(int argc, char *argv[])
{
  /* Different variables used for the optionparser as well as multiseq variable
     multiseq */
  MultiseqOptions options;

  try
  {
    options.parse(argc, argv);
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
  GttlMultiseq *multiseq = nullptr;
  RunTimeClass rt_multiseq{}, rt_total{};
  const std::vector<std::string> &inputfiles = options.inputfiles_get();
  try
  {
    const bool store_sequences = (options.width_option_get() >= 0 ||
                                  options.rankdist_option_is_set() ||
                                  options.short_header_option_is_set())
                                  ? true : false;
    if (options.zipped_option_is_set() && store_sequences)
    {
      multiseq = new GttlMultiseq(inputfiles[0],inputfiles[1],
                                  store_sequences,UINT8_MAX);
    } else
    {
      multiseq = new GttlMultiseq(inputfiles,store_sequences,UINT8_MAX);
    }
  }
  catch (std::string &msg)
  {
    for (auto &&inputfile : inputfiles)
    {
      std::cerr << argv[0] << ": file \"" << inputfile << "\""
                << msg << std::endl;
    }
    delete multiseq;
    return EXIT_FAILURE;
  }
  if (options.short_header_option_is_set())
  {
    multiseq->short_header_cache_create();
  }
  rt_multiseq.show("create GttlMultiseq");
  bool has_err = false;
  if (options.width_option_get() >= 0)
  {
    if (options.sample_size_get() > 0)
    {
      if (options.sample_size_get() > multiseq->sequences_number_get())
      {
        std::cerr << argv[0] << " you cannot sample "
                  << options.sample_size_get()
                  << " sequences from a set of "
                  << multiseq->sequences_number_get()
                  << " elements" << std::endl;
        has_err = true;
      } else
      {
        RunTimeClass rt_sample{};
        std::vector<size_t> sample
          = gttl_random_sample<size_t>(multiseq->sequences_number_get(),
                                       options.sample_size_get(),
                                       options.seed_get());
        for (auto &seqnum : sample)
        {
          multiseq->show_single_sequence(
                      static_cast<size_t>(options.width_option_get()),
                      options.short_header_option_is_set(),
                      seqnum);
        }
        StrFormat msg("create a sample of %lu (%.2f%%) from %lu sequences",
                      options.sample_size_get(),
                      100.0 * options.sample_size_get()/
                              multiseq->sequences_number_get(),
                      multiseq->sequences_number_get());
        rt_sample.show(msg.str());
      }
    } else
    {
      multiseq->show(static_cast<size_t>(options.width_option_get()),
                     options.short_header_option_is_set());
    }
  }
  if (options.statistics_option_is_set())
  {
    for (auto &&inputfile : inputfiles)
    {
      std::cout << "# filename\t" << inputfile << std::endl;
    }
    for (auto &msg : multiseq->statistics())
    {
      std::cout << "# " << msg << std::endl;
    }
  }
  if (options.rankdist_option_is_set())
  {
    if (options.protein_option_is_set())
    {
      static constexpr const char amino_acids[]
        = "A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y";
      LiterateMultiseq<amino_acids,20> literate_multiseq(multiseq);
      literate_multiseq.show_rank_dist(nullptr); /* no mutex necessary */
    } else
    {
      static constexpr const char nucleotides_upper_lower[] = "Aa|Cc|Gg|TtUu";
      LiterateMultiseq<nucleotides_upper_lower,4> literate_multiseq(multiseq);
      literate_multiseq.show_rank_dist(nullptr); /* no mutex necessary */
    }
  }
  delete multiseq;
  rt_total.show("total");
  return has_err ? EXIT_FAILURE : EXIT_SUCCESS;
}
