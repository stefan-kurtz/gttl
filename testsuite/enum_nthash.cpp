/*
  Copyright (c) 2021 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2021 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/
#include <cassert>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <string>
#include <algorithm>
#include <iostream>
#include "utilities/str_format.hpp"
#include "utilities/mathsupport.hpp"
#include "sequences/qgrams_hash_nthash.hpp"
#include "sequences/gttl_seq_iterator.hpp"
#include "sequences/guess_if_protein_seq.hpp"
#include "sequences/non_wildcard_ranges.hpp"

#ifndef NDEBUG
static void qgrams_nt_fwd_compare(alphabet::GttlAlphabet_UL_4 &alphabet,
                                  uint8_t *qgram_buffer,
                                  const char *orig_qgram,
                                  size_t qgram_length,
                                  uint64_t expected_hash_value)
{
  for (size_t idx = 0; idx < qgram_length; idx++)
  {
    qgram_buffer[idx] = alphabet.char_to_rank(orig_qgram[idx]);
  }
  uint64_t bf_hash_value = NTF64_first_hash_value_get(qgram_buffer,
                                                      qgram_length);
  assert(bf_hash_value == expected_hash_value);
}
#endif

static std::pair<uint64_t,size_t> apply_qgram_iterator(size_t qgram_length,
                                     const char *sequence,
                                     size_t this_length)
{
  uint8_t *qgram_buffer = new uint8_t [qgram_length];
#ifndef NDEBUG
  alphabet::GttlAlphabet_UL_4 dna_alphabet;
#endif
  uint64_t sum_hash_values = 0;
  size_t pos = 0;
  QgramNtHashFwdIterator4 qgiter(qgram_length,sequence,this_length);
  for (auto const &&code_pair : qgiter)
  {
    if (std::get<1>(code_pair) == 0)
    {
      sum_hash_values += std::get<0>(code_pair);
#ifndef NDEBUG
      qgrams_nt_fwd_compare(dna_alphabet,qgram_buffer,
                            sequence + pos,qgram_length,
                            std::get<0>(code_pair));
#endif
    }
    pos++;
  }
  delete [] qgram_buffer;
  return {sum_hash_values,pos};
}

static void enumerate_nt_hash_fwd(const char *inputfilename,size_t qgram_length)
{
  constexpr const int buf_size = 1 << 14;
  const bool is_protein = guess_if_protein_file(inputfilename);

  if (is_protein)
  {
    throw std::string(": can only handle DNA sequences");
    /* check_err.py checked */
  }
  GttlFpType in_fp = gttl_fp_type_open(inputfilename,"rb");
  if (in_fp == nullptr)
  {
    throw std::string(": cannot open file");
    /* check_err.py checked */
  }
  GttlSeqIterator<buf_size> gttl_si(in_fp);
  size_t total_length = 0;
  size_t num_of_sequences = 0;
  size_t max_sequence_length = 0;
  uint64_t sum_hash_values = 0;
  size_t count_all_qgrams = 0;
  try /* need this, as the catch needs to close the file pointer
         to prevent a memory leak */
  {
    for (auto &&si : gttl_si)
    {
      auto sequence = std::get<1>(si);
      total_length += sequence.size();
      max_sequence_length = std::max(max_sequence_length,sequence.size());
      std::vector<std::pair<size_t,size_t>> sequence_ranges;
      NonWildCardRangeIterator<'N'> nwcr_it(sequence.data(),sequence.size());
      sequence_ranges = nwcr_it.enumerate();
      for (auto &&nwr_it : sequence_ranges)
      {
        std::pair<uint64_t,size_t> result;
        const size_t this_length = std::get<1>(nwr_it) -
                                   std::get<0>(nwr_it) + 1;
        const char *seqptr = sequence.data() + std::get<0>(nwr_it);
        result = apply_qgram_iterator(qgram_length,seqptr,this_length);
        sum_hash_values += std::get<0>(result);
        count_all_qgrams += std::get<1>(result);
      }
      num_of_sequences++;
    }
  }
  catch (std::string &msg)
  {
    gttl_fp_type_close(in_fp);
    throw msg;
  }
  printf("# num_of_sequences\t%lu\n",num_of_sequences);
  printf("# total_length\t%lu\n",total_length);
  printf("# max_sequence_length\t%lu\n",max_sequence_length);
  printf("# num_of_sequences.bits\t%d\n",gt_required_bits(num_of_sequences));
  printf("# max_sequence_length.bits\t%d\n",
          gt_required_bits<size_t>(max_sequence_length));
  printf("# count_all_qgrams\t%lu\n",count_all_qgrams);
  printf("# sum_hash_values\t%lu\n",(unsigned long) sum_hash_values);
  gttl_fp_type_close(in_fp);
}

int main(int argc,char *argv[])
{
  const size_t qgram_length = 18;
  const char *progname = argv[0];
  bool haserr = false;
  for (int idx = 1; idx < argc; idx++)
  {
    try
    {
      enumerate_nt_hash_fwd(argv[idx],qgram_length);
    }
    catch (std::string &msg)
    {
      std::cerr << progname << ": file \"" << argv[idx] << "\""
                << msg << std::endl;
      haserr = true;
      break;
    }
  }
  return haserr ? EXIT_FAILURE : EXIT_SUCCESS;
}