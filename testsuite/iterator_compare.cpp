#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include "sequences/gttl_fasta_generator.hpp"
#include "sequences/gttl_fastq_generator.hpp"
#include "sequences/gttl_fastq_iterator.hpp"
#include "sequences/gttl_seq_iterator.hpp"
#include "utilities/cxxopts.hpp"
#include "utilities/gttl_line_iterator.hpp"
#include "utilities/runtime_class.hpp"

//TODO: valgrind
//TODO: Test with very-large-lines

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

  RunTimeClass runtime;


  if(result["fastq"].as<bool>())
  {
    GttlLineIterator<buf_size> line_it(file);
    GttlFastQIterator<GttlLineIterator<buf_size>> fq_it(line_it);

    runtime.reset();
    for(auto& entry : fq_it)
    {
      pseudo_hash_it += count_occ(entry.header_get().data(), '@');
      pseudo_hash_it += count_occ(entry.sequence_get().data(), 'a');
      pseudo_hash_it += count_occ(entry.quality_get().data(), 'F');
    }
    runtime.show("Iterator method");

    GttlFastQGenerator<buf_size> fq_gen(file);

    runtime.reset();
    for(auto entry : fq_gen)
    {
      pseudo_hash_gen += count_occ(entry->header, '@');
      pseudo_hash_gen += count_occ(entry->sequence, 'a');
      pseudo_hash_gen += count_occ(entry->quality, 'F');
    }
    runtime.show("Generator method");

    std::cout << pseudo_hash_it << "\t" << pseudo_hash_gen << "\n";
    // If the hashes are not identical, a mistake exists in
    // at least one implementation.
    // Performance comparinsons make no sense when this is the case.
    assert(pseudo_hash_it == pseudo_hash_gen);
  }
  else if(result["fasta"].as<bool>())
  {
    GttlSeqIterator<buf_size> fa_it(file);

    runtime.reset();
    for(auto &entry : fa_it)
    {
      pseudo_hash_it += count_occ(entry.header_get().data(), '>');
      pseudo_hash_it += count_occ(entry.sequence_get().data(), 'a');
      pseudo_hash_it += count_occ(entry.sequence_get().data(), 'A');
    }
    runtime.show("Iterator method");

    GttlFastAGenerator<buf_size> fa_gen(file);

    runtime.reset();
    for(auto entry : fa_gen)
    {
      pseudo_hash_gen += count_occ(entry->header, '>');
      pseudo_hash_gen += count_occ(entry->sequence, 'a');
      pseudo_hash_gen += count_occ(entry->sequence, 'A');
    }
    runtime.show("Generator method");

    std::cout << pseudo_hash_it << "\t" << pseudo_hash_gen << "\n";
    assert(pseudo_hash_it == pseudo_hash_gen);
  }
}
