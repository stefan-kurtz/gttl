#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstring>
#include <iostream>
#include "sequences/gttl_fasta_generator.hpp"
#include "sequences/gttl_fastq_generator.hpp"
#include "sequences/gttl_fastq_iterator.hpp"
#include "sequences/gttl_seq_iterator.hpp"
#include "utilities/cxxopts.hpp"
#include "utilities/gttl_line_iterator.hpp"

//TODO: Use utilities/RunTimeClass
//TODO: valgrind
//TODO: Test with very-large-lines

using hrc = std::chrono::high_resolution_clock;
using sec = std::chrono::seconds;

inline size_t count_occ(const char* string, char symbol)
{
  return std::count(string, string+std::strlen(string), symbol);
}

int main(int argc, char *argv[])
{
  cxxopts::Options options(argv[0],
                "Compare old Iterator vs new Generator implementations");
  options.set_width(80);
  options.set_tab_expansion();

  options.add_options()
      ("h,help", "Print usage information")
      ("f,file", "Input file to compare",
       cxxopts::value<std::string>())
      ("a,fasta", "Read FastA input")
      ("q,fastq", "Read FastQ input");


  options.parse_positional({"file"});
  options.positional_help("<input file>");

  auto result = options.parse(argc, argv);

  if(result.count("help"))
  {
    std::cout << options.help() << std::endl;
    exit(EXIT_SUCCESS);
  }

  if(result["fasta"].as<bool>() == result["fastq"].as<bool>())
  {
    std::cerr << argv[0] <<
        ": Exactly one of --fasta and --fastq must be passed.\n";
    exit(EXIT_FAILURE);
  }

  if(!result.count("file"))
  {
    std::cerr << argv[0] << ": A filename must be specified."
        << std::endl;
    exit(EXIT_FAILURE);
  }


  constexpr const size_t buf_size = (1 << 14);
  size_t pseudo_hash_gen = 0;
  size_t pseudo_hash_it = 0;

  const char* file = result["file"].as<std::string>().c_str();



  if(result["fastq"].as<bool>())
  {
    GttlLineIterator<buf_size> line_it(file);
    GttlFastQIterator<GttlLineIterator<buf_size>> fq_it(line_it);

    auto start_it = hrc::now();
    for(auto& entry : fq_it)
    {
      pseudo_hash_it += count_occ(entry.header_get().data(), '@');
      pseudo_hash_it += count_occ(entry.sequence_get().data(), 'a');
      pseudo_hash_it += count_occ(entry.quality_get().data(), 'F');
    }
    auto end_it = hrc::now();
    auto duration_it = std::chrono::duration_cast<sec>(end_it - start_it);
    std::cout << "Iterator method ran in:\t\t" <<
        duration_it.count() << " seconds.\n";

    GttlFastQGenerator<buf_size> fq_gen(file);

    auto start_gen = hrc::now();
    for(auto entry : fq_gen)
    {
      pseudo_hash_gen += count_occ(entry->header, '@');
      pseudo_hash_gen += count_occ(entry->sequence, 'a');
      pseudo_hash_gen += count_occ(entry->quality, 'F');
    }
    auto end_gen = hrc::now();
    auto duration_gen = std::chrono::duration_cast<sec>(end_gen-start_gen);

    // If the hashes are not identical, a mistake exists in
    // at least one implementation.
    // Performance comparinsons make no sense when this is the case.
    assert(pseudo_hash_it == pseudo_hash_gen);
    std::cout << "Generator method ran in:\t" <<
        duration_gen.count() << " seconds.\n";
  }
  else if(result["fasta"].as<bool>())
  {
    GttlSeqIterator<buf_size> fa_it(file);

    auto start_it = hrc::now();
    for(auto &entry : fa_it)
    {
      pseudo_hash_it += count_occ(entry.header_get().data(), '>');
      pseudo_hash_it += count_occ(entry.sequence_get().data(), 'a');
      pseudo_hash_it += count_occ(entry.sequence_get().data(), 'A');
    }
    auto end_it = hrc::now();
    auto duration_it = std::chrono::duration_cast<sec>(end_it - start_it);
    std::cout << "Iterator method ran in:\t\t" <<
      duration_it.count() << " seconds.\n";

    GttlFastAGenerator<buf_size> fa_gen(file);

    auto start_gen = hrc::now();
    for(auto entry : fa_gen)
    {
      pseudo_hash_gen += count_occ(entry->header, '>');
      pseudo_hash_gen += count_occ(entry->sequence, 'a');
      pseudo_hash_gen += count_occ(entry->sequence, 'A');
    }
    auto end_gen = hrc::now();
    auto duration_gen = std::chrono::duration_cast<sec>(end_gen - start_gen);

    assert(pseudo_hash_it == pseudo_hash_gen);
    std::cout << pseudo_hash_it << "\t" << pseudo_hash_gen << "\n";
    std::cout << "Generator method ran in:\t" <<
      duration_gen.count() << " seconds.\n";
  }
}
