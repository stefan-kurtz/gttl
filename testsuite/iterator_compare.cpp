#include <algorithm>
#include <cstring>
#include <iostream>
#include "sequences/gttl_fasta_generator.hpp"
#include "sequences/gttl_fastq_generator.hpp"
#include "sequences/gttl_fastq_iterator.hpp"
#include "sequences/gttl_seq_iterator.hpp"
#include "utilities/cxxopts.hpp"
#include "utilities/gttl_line_iterator.hpp"
#include "utilities/runtime_class.hpp"

inline static size_t count_occ(const std::string_view string, char symbol)
{
  return std::ranges::count(string, symbol);
}

inline static void assert_always(bool condition)
{
  if (not condition)
  {
    std::cerr <<
      "Assertion error: Iterators do not yield the same result!" << '\n';
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char *argv[])
{
  cxxopts::Options options(argv[0],"Compare old Iterator vs new Generator "
                                   "implementations");
  options.set_width(80);
  options.set_tab_expansion();

  options.add_options()
      ("h,help", "Print usage information")
      ("f,file", "Input file to compare",
       cxxopts::value<std::string>())
      ("g,generator_first", "first run generator")
      ("a,fasta", "Read FastA input")
      ("q,fastq", "Read FastQ input");

  options.parse_positional({"file"});
  options.positional_help("<input file>");

  auto result = options.parse(argc, argv);

  if (result.count("help"))
  {
    std::cout << options.help() << std::endl;
    exit(EXIT_SUCCESS);
  }

  if (result["fasta"].as<bool>() == result["fastq"].as<bool>())
  {
    std::cerr << argv[0] <<
        ": Exactly one of --fasta and --fastq must be passed.\n";
    exit(EXIT_FAILURE);
  }

  if (not result.count("file"))
  {
    std::cerr << argv[0] << ": A filename must be specified."
        << std::endl;
    exit(EXIT_FAILURE);
  }

  constexpr const size_t buf_size = (1 << 14);
  size_t pseudo_hash[3] = {0,0,0};

  const char* file = result["file"].as<std::string>().c_str();
  std::cout << "# " << file << '\n';

  RunTimeClass runtime{};
  const int first_idx = result["generator_first"].as<bool>() ? 1 : 0;
  if (result["fastq"].as<bool>())
  {
    for (int idx = 0; idx < 2; idx++)
    {
      if (idx == first_idx)
      {
        GttlLineIterator<buf_size> line_it(file);
        GttlFastQIterator<GttlLineIterator<buf_size>> fq_it(line_it);

        runtime.reset();
        for(auto& entry : fq_it)
        {
          pseudo_hash[1] += count_occ(entry.header_get(), '@');
          pseudo_hash[1] += count_occ(entry.sequence_get(), 'a');
          pseudo_hash[1] += count_occ(entry.quality_get(), 'F');
        }
        runtime.show("Iterator method");
      } else
      {
        GttlFastQGenerator<buf_size> fq_gen(file);
        runtime.reset();
        for(const auto *entry : fq_gen)
        {
          pseudo_hash[0] += count_occ(entry->header, '@');
          pseudo_hash[0] += count_occ(entry->sequence, 'a');
          pseudo_hash[0] += count_occ(entry->quality, 'F');
        }
        runtime.show("Generator method");
      }
    }
  } else
  {
    if (result["fasta"].as<bool>())
    {
      for (int idx = 0; idx < 2; idx++)
      {
        if (idx == first_idx)
        {
          GttlSeqIterator<buf_size> fa_it(file);
          runtime.reset();
          for(auto &entry : fa_it)
          {
            pseudo_hash[idx] += count_occ(entry.header_get(), ':');
            pseudo_hash[idx] += count_occ(entry.sequence_get(), 'a');
            pseudo_hash[idx] += count_occ(entry.sequence_get(), 'A');
          }
          runtime.show("Iterator method");
        } else
        {
          GttlFastAGenerator<buf_size> fa_gen(file);
          runtime.reset();
          for(const auto *entry : fa_gen)
          {
            pseudo_hash[idx] += count_occ(entry->header_get(), ':');
            pseudo_hash[idx] += count_occ(entry->sequence_get(), 'a');
            pseudo_hash[idx] += count_occ(entry->sequence_get(), 'A');
          }
          runtime.show("Generator method");
        }
      }
    }
  }
  for (int idx = 0; idx < 2; idx++)
  {
    std::cout << pseudo_hash[idx] << (idx == 0 ? "\t" : "\n");
  }
  // If the hashes are not identical, a mistake exists in
  // at least one implementation.
  // Performance comparinsons make no sense when this is the case.
  assert_always(pseudo_hash[0] == pseudo_hash[1]);
  return EXIT_SUCCESS;
}
