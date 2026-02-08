#include <format>
#include <iostream>
#include <cstdint>
#include <cinttypes>
#include <vector>
#include "utilities/runtime_class.hpp"
#include "utilities/gttl_line_generator.hpp"
#include "sequences/gttl_fasta_generator.hpp"
#include "sequences/myersapm.hpp"

/* for all patterns for all sequences */
static void myers_bitvector_algorithm_p_times_s(
  const std::string &pattern_file,
  const std::string &sequence_file,
  size_t error_threshold)
{
  RunTimeClass total{};
  GttlLineGenerator gttl_lg(pattern_file);
  size_t count_match_events = 0;
  size_t count_patterns = 0;
  for (const auto &pattern : gttl_lg)
  {
    if (error_threshold > pattern.size()/2)
    {
      throw std::format("error threshold {} is too large; maximum value for "
                        "pattern of length {} is {}", error_threshold,
                        pattern.size(), pattern.size()/2);
    }
    count_patterns++;
    MyersBitvectorAlgorithm<uint64_t> myers_bitvector_algorithm(pattern);
    constexpr const size_t buf_size = size_t{1} << size_t{14};
    GttlFastAGenerator<buf_size> gttl_si(sequence_file.c_str());
    size_t count_matches = 0;
    for (auto &&si : gttl_si)
    {
      count_match_events++;
      const std::string_view &sequence = si->sequence_get();
      for (size_t j = 0; j < sequence.size(); j++)
      {
        const size_t cost = myers_bitvector_algorithm.transform(sequence[j]);
        if (cost <= error_threshold)
        {
          count_matches++;
        }
      }
      myers_bitvector_algorithm.reset_for_next_sequence();
    }
    std::cout << pattern << '\t' << count_matches << '\n';
  }
  assert(count_match_events % count_patterns == 0);
  total.show(std::format("p_times_s: match {} patterns against {} sequences",
                         count_patterns, count_match_events/count_patterns));
}

/* for all sequences for all patterns */
static void myers_bitvector_algorithm_s_times_p(
  const std::string &pattern_file,
  const std::string &sequence_file,
  size_t error_threshold)
{
  RunTimeClass total{};
  GttlLineGenerator gttl_lg(pattern_file);
  std::vector<std::string> pattern_vector;
  std::vector<MyersBitvectorAlgorithm<uint64_t>> mbv_vector;
  for (const auto &pattern : gttl_lg)
  {
    if (error_threshold > pattern.size()/2)
    {
      throw std::format("error threshold {} is too large; maximum value for "
                        "pattern of length {} is {}", error_threshold,
                        pattern.size(), pattern.size()/2);
    }
    pattern_vector.push_back(pattern);
    mbv_vector.push_back(MyersBitvectorAlgorithm<uint64_t>(pattern));
  }
  constexpr const size_t buf_size = size_t{1} << size_t{14};
  GttlFastAGenerator<buf_size> gttl_si(sequence_file.c_str());
  std::vector<size_t> count_matches_vector(mbv_vector.size(), 0);
  size_t count_sequences = 0;
  for (auto &&si : gttl_si)
  {
    count_sequences++;
    const std::string_view &sequence = si->sequence_get();
    for (size_t p_idx = 0; p_idx < mbv_vector.size(); p_idx++)
    {
      MyersBitvectorAlgorithm<uint64_t> &myers_bitvector_algorithm
        = mbv_vector[p_idx];
      for (size_t j = 0; j < sequence.size(); j++)
      {
        const size_t cost = myers_bitvector_algorithm.transform(sequence[j]);
        if (cost <= error_threshold)
        {
          count_matches_vector[p_idx]++;
        }
      }
      myers_bitvector_algorithm.reset_for_next_sequence();
    }
  }
  for (size_t p_idx = 0; p_idx < count_matches_vector.size(); p_idx++)
  {
    std::cout << pattern_vector[p_idx] << '\t' << count_matches_vector[p_idx]
              << '\n';
  }
  total.show(std::format("s_times_p: match {} sequences against all {} "
                         "patterns", count_sequences, pattern_vector.size()));
}

static void usage(const char *progname)
{
  std::cerr << "Usage: " << progname
            << (" p_times_s|s_times_p <num_threads> <pattern_file> "
                "<sequence_file> <error_threshold>\n");
}

int main(int argc, char *argv[])
{
  if (argc != 6 or (std::string(argv[1]) != "p_times_s" and
                    std::string(argv[1]) != "s_times_p"))
  {
    usage(argv[0]);
    return EXIT_FAILURE;
  }
  int64_t read_int;
  if (sscanf(argv[5], "%" PRIi64, &read_int) != 1 or read_int < 0)
  {
    usage(argv[0]);
    std::cerr << "Cannot read non-negative integer from " << argv[5] << '\n';
    return EXIT_FAILURE;
  }
  const size_t error_threshold = static_cast<size_t>(read_int);
  if (sscanf(argv[2], "%" PRIi64, &read_int) != 1 or read_int <= 0)
  {
    usage(argv[0]);
    std::cerr << "Cannot read positive integer from " << argv[2] << '\n';
    return EXIT_FAILURE;
  }
  [[maybe_unused]] const size_t num_threads = static_cast<size_t>(read_int);
  const std::string pattern_file{argv[3]};
  const std::string sequence_file{argv[4]};
  try
  {
    (std::string(argv[1]) == "p_times_s"
      ? myers_bitvector_algorithm_p_times_s
      : myers_bitvector_algorithm_s_times_p) (pattern_file,
                                              sequence_file,
                                              error_threshold);
  }
  catch (const std::runtime_error &err)
  {
    std::cerr << argv[0] << ": " <<  err.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
