#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <cinttypes>
#include <stdexcept>
#include <string>
#include "utilities/gttl_file_open.hpp"
#include "sequences/char_range.hpp"
#include "sequences/char_finder.hpp"
#include "sequences/gttl_fasta_generator.hpp"
#include "sequences/qgrams_hash_invint.hpp"

static constexpr const char_finder::NucleotideFinder nucleotide_finder{};

#ifndef NDEBUG
template<class HashValuePairIterator>
static void verify_hash_value_pair(HashValuePairIterator &qgiter,
                                   uint64_t hash_value,
                                   uint64_t compl_hash_value)
{
  size_t qgram_length = qgiter.qgram_length_get();
  const uint8_t *qgram = qgiter.qgram_decode(hash_value);
  uint8_t *qgram_direct = new uint8_t [qgram_length];
  memcpy(qgram_direct,qgram,qgram_length * sizeof *qgram);
  const uint8_t *qgram_rc = qgiter.qgram_decode(compl_hash_value);
  for (size_t idx = 0; idx < qgram_length; idx++)
  {
    uint8_t cc = qgram_rc[qgram_length - 1 - idx];
    assert(cc < 4);
    cc = complement_uint8(cc);
    if (cc != qgram_direct[idx])
    {
      delete[] qgram_direct;
      throw std::runtime_error(std::string("incorrect reverse complement "
                                           "hash_value=") +
                               std::to_string(hash_value) +
                               std::string("\t") +
                               std::string("compl_hash_value=") +
                               std::to_string(compl_hash_value));
    }
  }
  delete[] qgram_direct;
}
#endif

template<class HashValuePairIterator>
static void verify_hashvalues_for_file(const char *inputfilename,
                                       size_t qgram_length)
{
  const GttlFpType in_fp = gttl_fp_type_open(inputfilename, "rb");
  using NucleotideRanger = GttlCharRange<char_finder::NucleotideFinder,
                                         nucleotide_finder,
                                         true,false>;
  if (in_fp == nullptr)
  {
    throw std::runtime_error(std::string("cannot open file ") + inputfilename);
  }
  GttlFastAGenerator fasta_gen(in_fp);
    size_t ranges_total_length = 0;
    size_t hash_value_sum = 0;
    size_t compl_hash_value_sum = 0;
  for (const auto *si : fasta_gen)
  {
    auto sequence = si->sequence_get();
    const NucleotideRanger nuc_ranger(sequence.data(), sequence.size());
    for (auto const &&range : nuc_ranger)
    {
      const size_t this_length = std::get<1>(range);
      const char *const substring = sequence.data() + std::get<0>(range);
      HashValuePairIterator qgiter(qgram_length,substring, this_length);
      for (auto const &&code_pair : qgiter)
      {
        const uint64_t hash_value       = std::get<0>(code_pair);
        const uint64_t compl_hash_value = std::get<1>(code_pair);
        hash_value_sum += hash_value;
        compl_hash_value_sum += compl_hash_value;
#ifndef NDEBUG
        verify_hash_value_pair<HashValuePairIterator>
                              (qgiter,hash_value,compl_hash_value);
#endif
      }
      ranges_total_length += this_length;
    }
  }
  std::cout << "# ranges_total_length\t" << ranges_total_length << '\n';
  std::cout << "# hash_value_sum\t" << hash_value_sum << '\n';
  std::cout << "# complement hash_value_sum\t" << compl_hash_value_sum << '\n';
}

int main(int argc,char *argv[])
{
  int64_t read_long;
  if (argc != 3 || sscanf(argv[1],"%" PRIi64,&read_long) != 1 || read_long <= 0)
  {
    std::cerr << "Usage: " << argv[0] << ": <kmer_length> <inputfile>\n";
    return EXIT_FAILURE;
  }
  const size_t qgram_length = static_cast<size_t>(read_long);
  const char *const inputfilename = argv[2];
  try
  {
    verify_hashvalues_for_file<InvertibleIntegercode2Iterator4>
                              (inputfilename,qgram_length);
  }
  catch (const std::exception &err)
  {
    std::cerr << argv[0] << ": " << err.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
