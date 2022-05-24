#include <cstddef>
#include <cstdio>
#include <iostream>
#include "sequences/char_range.hpp"
#include "sequences/char_finder.hpp"
#include "sequences/gttl_seq_iterator.hpp"
#include "sequences/qgrams_hash2_invint.hpp"

static constexpr const char_finder::NucleotideFinder nucleotide_finder{};

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
  using NucleotideRanger = GttlCharRange<char_finder::NucleotideFinder,
                                         nucleotide_finder,
                                         true,false>;
  try
  {
    constexpr const int buf_size = 1 << 14;
    GttlSeqIterator<buf_size> gttl_si(argv[2]);
    size_t ranges_total_length = 0, seqnum = 0;
    for (auto &&si : gttl_si)
    {
      auto sequence = si.sequence_get();
      NucleotideRanger nuc_ranger(sequence.data(),sequence.size());
      for (auto const &&range : nuc_ranger)
      {
        const size_t this_length = std::get<1>(range);
        const char *substring = sequence.data() + std::get<0>(range);
        std::cout << seqnum << "\t" << std::get<0>(range)
                  << "\t" << this_length << std::endl;
        InvertibleIntegercode2Iterator4 qgiter(qgram_length,substring,
                                               this_length);
        for (auto const &&code_pair : qgiter)
        {
          std::cout << std::get<0>(code_pair) << "\t" << std::get<1>(code_pair)
                    << std::endl;
        }
        ranges_total_length += this_length;
      }
      seqnum++;
    }
    std::cout << "# ranges_total_length\t" << ranges_total_length << std::endl;
  }
  catch (std::string &msg)
  {
    throw msg;
  }
  return EXIT_SUCCESS;
}
