#include <format>
#include <iostream>
#include <cstdint>
#include <cinttypes>
#include "utilities/gttl_line_generator.hpp"
#include "sequences/gttl_fasta_generator.hpp"
#include "sequences/myersapm.hpp"

static void run_myers_bitvector_algorithm(const std::string &pattern_file,
                                          const std::string &sequence_file,
                                          size_t error_threshold)
{
  GttlLineGenerator gttl_lg(pattern_file);
  for (const auto &pattern : gttl_lg)
  {
    if (error_threshold > pattern.size()/2)
    {
      throw std::format("error threshold {} is too large; maximum value for "
                        "pattern of length {} is {}", error_threshold,
                        pattern.size(), pattern.size()/2);
    }
    MyersBitvectorAlgorithm<uint64_t> myers_bitvector_algorithm(pattern);
    constexpr const size_t buf_size = size_t{1} << size_t{14};
    GttlFastAGenerator<buf_size> gttl_si(sequence_file.c_str());
    for (auto &&si : gttl_si)
    {
      const std::string_view &sequence = si->sequence_get();
      size_t count_matches = 0;
      for (size_t j = 0; j < sequence.size(); j++)
      {
        const size_t cost = myers_bitvector_algorithm.transform(sequence[j]);
        if (cost <= error_threshold)
        {
          count_matches++;
        }
      }
      std::cout << pattern << '\t' << count_matches << '\n';
      myers_bitvector_algorithm.reset_for_next_sequence();
    }
  }
}

int main(int argc,char *argv[])
{
  if (argc != 4)
  {
    std::cerr << "Usage: " << argv[0]
              << " <pattern_file> <sequence_file> <error_threshold>\n";
    return EXIT_FAILURE;
  }
  int64_t read_int;
  if (sscanf(argv[3],"%" PRIi64,&read_int) != 1 or read_int < 0)
  {
    std::cerr << "Usage: " << argv[0]
              << " <pattern_file> <sequence_file> <error_threshold>\n";
    std::cerr << "Cannot read non-negative integer from " << argv[3] << '\n';
    return EXIT_FAILURE;
  }
  const size_t error_threshold = static_cast<size_t>(read_int);
  const std::string pattern_file{argv[1]};
  const std::string sequence_file{argv[2]};
  try
  {
    run_myers_bitvector_algorithm(pattern_file,
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
