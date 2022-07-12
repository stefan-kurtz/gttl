#include <cstddef>
#include <cstdio>
#include <iostream>
#include "sequences/char_range.hpp"
#include "sequences/char_finder.hpp"
#include "sequences/gttl_seq_iterator.hpp"
#include "sequences/qgrams_hash_invint.hpp"

static constexpr const char_finder::NucleotideFinder nucleotide_finder{};

#ifndef NDEBUG
static void verify_hash_values(const char *progname,
                               InvertibleIntegercode2Iterator4 &qgiter,
                               uint64_t hash_value,
                               uint64_t compl_hash_value)
{
  static constexpr const uint8_t complement[] = {uint8_t(2),uint8_t(3),
                                                 uint8_t(0),uint8_t(1)};
  size_t qgram_length = qgiter.qgram_length_get();
  const uint8_t *qgram = qgiter.qgram_decode(hash_value);
  uint8_t *qgram_direct = new uint8_t [qgram_length];
  memcpy(qgram_direct,qgram,qgram_length * sizeof *qgram);
  const uint8_t *qgram_rc = qgiter.qgram_decode(compl_hash_value);
  for (size_t idx = 0; idx < qgram_length; idx++)
  {
    uint8_t cc = qgram_rc[qgram_length - 1 - idx];
    assert(cc < 4);
    cc = complement[cc];
    if (cc != qgram_direct[idx])
    {
      std::cerr << progname << ": incorrect reverse complement integer code"
                << "hash_value=" << hash_value
                << "\tcompl_hash_value=" << compl_hash_value << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  delete[] qgram_direct;
}
#endif

int main(int argc,char *argv[])
{
  long read_long;
  if (argc != 3 || sscanf(argv[1],"%ld",&read_long) != 1 || read_long <= 0)
  {
    std::cerr << "Usage: " << argv[0] << ": <kmer_length> <inputfile>"
              << std::endl;
    return EXIT_FAILURE;
  }
  const size_t qgram_length = static_cast<size_t>(read_long);
  const char *inputfilename = argv[2];
  using NucleotideRanger = GttlCharRange<char_finder::NucleotideFinder,
                                         nucleotide_finder,
                                         true,false>;
  constexpr const int buf_size = 1 << 14;
  GttlFpType in_fp = gttl_fp_type_open(inputfilename,"rb");
  if (in_fp == nullptr)
  {
    std::cerr << argv[0] << " cannot open file " << inputfilename << std::endl;
    return EXIT_FAILURE;
  }
  GttlSeqIterator<buf_size> gttl_si(in_fp);
  try /* need this, as the catch needs to close the file pointer
         to prevent a memory leak */
  {
    size_t ranges_total_length = 0,
           hash_value_sum = 0,
           compl_hash_value_sum = 0;
#ifndef NDEBUG
    size_t seqnum = 0;
#endif
    for (auto &&si : gttl_si)
    {
      auto sequence = si.sequence_get();
      NucleotideRanger nuc_ranger(sequence.data(),sequence.size());
      for (auto const &&range : nuc_ranger)
      {
        const size_t this_length = std::get<1>(range);
        const char *substring = sequence.data() + std::get<0>(range);
        InvertibleIntegercode2Iterator4 qgiter(qgram_length,substring,
                                               this_length);
        for (auto const &&code_pair : qgiter)
        {
          uint64_t hash_value = std::get<0>(code_pair),
                   compl_hash_value = std::get<1>(code_pair);
          hash_value_sum += hash_value;
          compl_hash_value_sum += compl_hash_value;
#ifndef NDEBUG
          verify_hash_values(argv[0],qgiter,hash_value,compl_hash_value);
#endif
        }
        ranges_total_length += this_length;
      }
#ifndef NDEBUG
      seqnum++;
#endif
    }
    std::cout << "# ranges_total_length\t" << ranges_total_length << std::endl;
    std::cout << "# hash_value_sum\t" << hash_value_sum << std::endl;
    std::cout << "# complement hash_value_sum\t" << compl_hash_value_sum
              << std::endl;
  }
  catch (std::string &msg)
  {
    gttl_fp_type_close(in_fp);
    std::cerr << argv[0] << msg << std::endl;
    return EXIT_FAILURE;
  }
  gttl_fp_type_close(in_fp);
  return EXIT_SUCCESS;
}
