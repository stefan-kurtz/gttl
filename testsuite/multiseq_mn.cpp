#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <limits>
#include <string>
#include <vector>
#include <format>
#include "sequences/gttl_multiseq.hpp"
#include "utilities/cxxopts.hpp"
#include "utilities/runtime_class.hpp"
#include "sequences/literate_multiseq.hpp"
#include "utilities/random_sample.hpp"

static void usage(const cxxopts::Options &options)
{
  std::cerr << options.help() << '\n';
}

class MultiseqOptions
{
 private:
  std::vector<std::string> inputfiles;
  bool help_option;
  bool protein_option;
  bool zipped_option;
  bool rankdist_option;
  bool short_header_option;
  bool sorted_by_header_option;
  bool statistics_option;
  size_t sample_size;
  size_t length_dist_bin_size;
  size_t min_length;
  size_t max_length;
  unsigned int seed;
  int width_arg = -1;
 public:
  MultiseqOptions(void)
   : help_option(false)
   , protein_option(false)
   , zipped_option(false)
   , rankdist_option(false)
   , short_header_option(false)
   , sorted_by_header_option(false)
   , statistics_option(false)
   , sample_size(0)
   , length_dist_bin_size(0)
   , min_length(0)
   , max_length(0)
   , seed(0)
 { }

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
       ("l,length_dist", "output distribution of lengths of sequences; "
                         "argument specifies the size of the bins",
        cxxopts::value<size_t>(length_dist_bin_size)->default_value("0"))
       ("min_length", "specify minimum length of sequences to select when "
                      "sampling or formatting sequences",
        cxxopts::value<size_t>(min_length)->default_value("0"))
       ("max_length", "specify maximum length of sequences to select when "
                      "sampling or formatting sequences",
        cxxopts::value<size_t>(max_length)->default_value("0"))
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
       ("sorted_by_header", "output sequences lexicographically sorted by "
                            "the header; if option -s/--short_header is used, "
                            "then only the short header determines the order; "
                            "option requires to use option -w/--width",
        cxxopts::value<bool>(sorted_by_header_option)->default_value("false"))
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
        for (const auto & unmatched_arg : unmatched_args)
        {
          inputfiles.push_back(unmatched_arg);
        }
      }
      if (zipped_option && inputfiles.size() != 2)
      {
        throw std::invalid_argument("option -z/--zipped requires exactly "
                                    "two files");
      }
      if (sorted_by_header_option and width_arg == -1)
      {
        throw std::invalid_argument("option --sorted_by_header requires to "
                                    "use option -w/--width");
      }
      if (length_dist_bin_size > 0)
      {
        if (sample_size > 0)
        {
          throw std::invalid_argument("options --sample and --length_dist "
                                      "are not compatible");
        }
        if (seed > 0)
        {
          throw std::invalid_argument("options --seed and --length_dist "
                                      "are not compatible");
        }
        if (zipped_option)
        {
          throw std::invalid_argument("options --zipped and --length_dist "
                                      "are not compatible");
        }
        if (width_arg > 0)
        {
          throw std::invalid_argument("options --width and --length_dist "
                                      "are not compatible");
        }
        if (min_length > 0)
        {
          throw std::invalid_argument("options --min_length and --length_dist "
                                      "are not compatible");
        }
        if (max_length > 0)
        {
          throw std::invalid_argument("options --max_length and --length_dist "
                                      "are not compatible");
        }
        if (short_header_option)
        {
          throw std::invalid_argument("options --short_header and "
                                      "--length_dist are not compatible");
        }
      }
    }
    catch (const cxxopts::exceptions::exception &e)
    {
      usage(options);
      throw std::invalid_argument(e.what());
    }
  }
  [[nodiscard]] bool help_option_is_set(void) const noexcept
  {
    return help_option;
  }
  [[nodiscard]] bool protein_option_is_set(void) const noexcept
  {
    return protein_option;
  }
  [[nodiscard]] bool zipped_option_is_set(void) const noexcept
  {
    return zipped_option;
  }
  [[nodiscard]] bool rankdist_option_is_set(void) const noexcept
  {
    return rankdist_option;
  }
  [[nodiscard]] bool short_header_option_is_set(void) const noexcept
  {
    return short_header_option;
  }
  [[nodiscard]] bool sorted_by_header_option_is_set(void) const noexcept
  {
    return sorted_by_header_option;
  }
  [[nodiscard]] bool statistics_option_is_set(void) const noexcept
  {
    return statistics_option;
  }
  [[nodiscard]] int width_option_get(void) const noexcept
  {
    return width_arg;
  }
  [[nodiscard]] const std::vector<std::string> &inputfiles_get(void) const
    noexcept
  {
    return inputfiles;
  }
  [[nodiscard]] size_t sample_size_get(void) const noexcept
  {
    return sample_size;
  }
  [[nodiscard]] size_t min_length_get(void) const noexcept
  {
    return min_length;
  }
  [[nodiscard]] size_t max_length_get(void) const noexcept
  {
    return max_length;
  }
  [[nodiscard]] unsigned int seed_get(void) const noexcept
  {
    return seed;
  }
  [[nodiscard]] size_t length_dist_bin_size_get(void) const noexcept
  {
    return length_dist_bin_size;
  }
};

static std::vector<size_t> min_max_length_random_sample(
                                                   size_t *min_max_satisfy,
                                                   const GttlMultiseq *multiseq,
                                                   size_t min_length,
                                                   size_t max_length,
                                                   size_t sample_size,
                                                   unsigned int seed)
{
  if (max_length == 0)
  {
    max_length = std::numeric_limits<size_t>::max();
  }
  std::vector<size_t> selected_seqnums;
  for (size_t seqnum = 0; seqnum < multiseq->sequences_number_get(); seqnum++)
  {
    const size_t current_length = multiseq->sequence_length_get(seqnum);
    if (current_length >= min_length and current_length <= max_length)
    {
      selected_seqnums.push_back(seqnum);
    }
  }
  *min_max_satisfy = selected_seqnums.size();
  const std::vector<size_t> index_sample
    = gttl_random_sample<size_t>(selected_seqnums.size(),
                                 sample_size,
                                 seed);
  std::vector<size_t> sample;
  sample.reserve(sample_size);
  for (const auto idx : index_sample)
  {
    sample.push_back(selected_seqnums[idx]);
  }
  return sample;
}

int main(int argc, char *argv[])
{
  /* Different variables used for the optionparser as well as multiseq variable
     multiseq */
  MultiseqOptions options;

  try
  {
    options.parse(argc, argv);
  }
  catch (std::invalid_argument &err) /* check_err.py */
  {
    std::cerr << argv[0] << ": " << err.what() << '\n';
    return EXIT_FAILURE;
  }
  if (options.help_option_is_set())
  {
    return EXIT_SUCCESS;
  }
  GttlMultiseq *multiseq = nullptr;
  RunTimeClass rt_multiseq{};
  RunTimeClass rt_total{};
  const std::vector<std::string> &inputfiles = options.inputfiles_get();
  try
  {
    constexpr const bool store_header = true;
    const bool store_sequence = (options.width_option_get() >= 0 ||
                                 options.rankdist_option_is_set() ||
                                 options.short_header_option_is_set())
                                 ? true : false;
    const uint8_t padding_char = UINT8_MAX;
    if (options.zipped_option_is_set() && store_sequence)
    {
      multiseq = new GttlMultiseq(inputfiles[0],inputfiles[1],
                                  store_header,store_sequence,padding_char);
    } else
    {
      multiseq = new GttlMultiseq(inputfiles,store_header,store_sequence,
                                  padding_char);
    }
  }
  catch (const std::exception &err)
  {
    for (auto &&inputfile : inputfiles)
    {
      std::cerr << argv[0] << ": file \"" << inputfile << "\"" << err.what()
                << '\n';
    }
    delete multiseq;
    return EXIT_FAILURE;
  }
  if (options.short_header_option_is_set())
  {
    multiseq->short_header_cache_create();
  }
  rt_multiseq.show("create GttlMultiseq");

  if (options.width_option_get() >= 0)
  {
    if (options.sample_size_get() > 0)
    {
      try
      {
        size_t min_max_satisfy = 0;
        RunTimeClass rt_sample{};
        const std::vector<size_t> sample
          = (options.min_length_get() > 0 or options.max_length_get() > 0)
              ? min_max_length_random_sample(&min_max_satisfy,
                                             multiseq,
                                             options.min_length_get(),
                                             options.max_length_get(),
                                             options.sample_size_get(),
                                             options.seed_get())
              : gttl_random_sample<size_t>(multiseq->sequences_number_get(),
                                           options.sample_size_get(),
                                           options.seed_get());
        for (const auto seqnum : sample)
        {
          multiseq->show_single_sequence(
                      static_cast<size_t>(options.width_option_get()),
                      options.short_header_option_is_set(),
                      seqnum);
        }
        const std::string msg
          = (options.min_length_get() > 0 or options.max_length_get() > 0)
            ? std::format("create a sample of {} ({:.2f}%) "
                          "from {} sequences satisfying the min_length/"
                          "max_length constraints from a total number of "
                          "{} sequences",
                          options.sample_size_get(),
                          100.0 * options.sample_size_get()/min_max_satisfy,
                          min_max_satisfy,
                          multiseq->sequences_number_get())
            : std::format("create a sample of {} ({:.2f}%%) from {} sequences",
                          options.sample_size_get(),
                          100.0 * options.sample_size_get()
                                  / multiseq->sequences_number_get(),
                          multiseq->sequences_number_get());
        rt_sample.show(msg);
      }
      catch (const std::runtime_error &msg)
      {
        std::cerr << argv[0] << ": " << msg.what() << '\n';
        return EXIT_FAILURE;
      }
    } else
    {
      if (options.min_length_get() > 0 or options.max_length_get() > 0)
      {
        const size_t min_length = options.min_length_get();
        const size_t max_length = options.max_length_get() == 0
                                    ? std::numeric_limits<size_t>::max()
                                    : options.max_length_get();
        for (size_t seqnum = 0; seqnum < multiseq->sequences_number_get();
             seqnum++)
        {
          const size_t current_length = multiseq->sequence_length_get(seqnum);
          if (current_length >= min_length and current_length <= max_length)
          {
            multiseq->show_single_sequence(
                        static_cast<size_t>(options.width_option_get()),
                        options.short_header_option_is_set(),
                        seqnum);
          }
        }
      } else
      {
        if (options.sorted_by_header_option_is_set())
        {
          multiseq->show_sorted_by_header(
                      static_cast<size_t>(options.width_option_get()),
                      options.short_header_option_is_set());
        } else
        {
          multiseq->show(static_cast<size_t>(options.width_option_get()),
                         options.short_header_option_is_set());
        }
      }
    }
  }
  if (options.statistics_option_is_set())
  {
    for (auto &&inputfile : inputfiles)
    {
      std::cout << "# filename\t" << inputfile << '\n';
    }
    for (auto &msg : multiseq->statistics())
    {
      std::cout << "# " << msg << '\n';
    }
  }
  if (options.rankdist_option_is_set())
  {
    if (options.protein_option_is_set())
    {
      static constexpr const char amino_acids[]
        = "A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y";
      const LiterateMultiseq<amino_acids, 20> literate_multiseq(multiseq);
      literate_multiseq.show_rank_dist(nullptr); /* no mutex necessary */
    } else
    {
      static constexpr const char nucleotides_upper_lower[] = "Aa|Cc|Gg|TtUu";
      const LiterateMultiseq<nucleotides_upper_lower, 4>
           literate_multiseq(multiseq);
      literate_multiseq.show_rank_dist(nullptr); /* no mutex necessary */
    }
  }
  if (options.length_dist_bin_size_get() > 0)
  {
    const size_t bin_size = options.length_dist_bin_size_get();
    const auto length_dist = multiseq->length_distribution(bin_size);
    for (const auto &[length, count] : length_dist)
    {
      std::cout << (bin_size * length) << "--" << (bin_size * (length + 1) - 1)
                << "\t" << count << '\n';
    }
  }
  delete multiseq;
  rt_total.show("total");
  return EXIT_SUCCESS;
}
