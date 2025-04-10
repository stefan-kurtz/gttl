#include <algorithm>
#include <cstring>
#include <iostream>
#include "sequences/gttl_fasta_generator.hpp"
#include "sequences/gttl_fastq_generator.hpp"
#include "sequences/gttl_fastq_iterator.hpp"
#include "sequences/gttl_seq_iterator.hpp"
#include "utilities/cxxopts.hpp"
#include "utilities/constexpr_for.hpp"
#include "utilities/gttl_line_iterator.hpp"
#include "utilities/runtime_class.hpp"

//TODO: Test with very-large-lines

inline static size_t count_occ(const std::string_view string, char symbol)
{
  return std::ranges::count(string, symbol);
}

inline static void assert_always(bool condition)
{
  if (not condition)
  {
    std::cerr <<
      "Assertion error: Iterators do not yield the same result!" << std::endl;
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

  if (!result.count("file"))
  {
    std::cerr << argv[0] << ": A filename must be specified."
        << std::endl;
    exit(EXIT_FAILURE);
  }

  constexpr const size_t buf_size = (1 << 14);
  size_t pseudo_hash[3] = {0,0,0};

  const char* file = result["file"].as<std::string>().c_str();

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
          pseudo_hash[idx] += count_occ(entry.header_get().data(), '@');
          pseudo_hash[idx] += count_occ(entry.sequence_get().data(), 'a');
          pseudo_hash[idx] += count_occ(entry.quality_get().data(), 'F');
        }
        runtime.show("Iterator method");
      } else
      {
        constexpr_for<1,2+1,1>([&](auto mode)
        {
          constexpr const bool use_heap = mode == 1 ? false : true;
          GttlFastQGenerator<buf_size, use_heap> fq_gen(file);
          runtime.reset();
          int count_idx;
          if (idx == 1)
          {
            count_idx = mode;
          } else
          {
            count_idx = (mode == 1) ? 0 : 2;
          }
          for(auto entry : fq_gen)
          {
            pseudo_hash[count_idx] += count_occ(entry->header, '@');
            pseudo_hash[count_idx] += count_occ(entry->sequence, 'a');
            pseudo_hash[count_idx] += count_occ(entry->quality, 'F');
          }
          runtime.show(use_heap ? "Generator method with heap" :
                                  "Generator method with stack");
        });
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
          constexpr const bool use_heap = true;
          GttlFastAGenerator<buf_size, use_heap> fa_gen(file);
          runtime.reset();
          for(auto entry : fa_gen)
          {
            pseudo_hash[idx] += count_occ(entry->header_get(), ':');
            pseudo_hash[idx] += count_occ(entry->sequence_get(), 'a');
            pseudo_hash[idx] += count_occ(entry->sequence_get(), 'A');
          }
          runtime.show("Generator method with heap");
        }
      }
    }
  }
  const int num_values = result["fasta"].as<bool>() ? 2 : 3;
  for (int idx = 0; idx < num_values; idx++)
  {
    std::cout << pseudo_hash[idx] << (idx < (num_values - 1) ? "\t" : "\n");
  }
  // If the hashes are not identical, a mistake exists in
  // at least one implementation.
  // Performance comparinsons make no sense when this is the case.
  assert_always(pseudo_hash[0] == pseudo_hash[1]);
  if (num_values > 2)
  {
    assert_always(pseudo_hash[1] == pseudo_hash[2]);
  }
  return EXIT_SUCCESS;
}
