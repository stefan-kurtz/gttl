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
#include <ios>
#include "utilities/mathsupport.hpp"
#include "sequences/qgrams_hash_nthash.hpp"
#include "sequences/gttl_seq_iterator.hpp"
#include "sequences/guess_if_protein_seq.hpp"
#include "sequences/non_wildcard_ranges.hpp"

#ifndef NDEBUG
template<char wildcard>
static void verify_non_wildcard_ranges(const std::string &sequence,
            NonWildCardRangeVector &non_wildcard_ranges)
{
  size_t wildcard_start = 0;
  for (auto &&nwr_it : non_wildcard_ranges)
  {
    for (size_t idx = wildcard_start; idx < std::get<0>(nwr_it); idx++)
    {
      assert(sequence.at(idx) == wildcard);
    }
    for (size_t idx = std::get<0>(nwr_it); idx <= std::get<1>(nwr_it); idx++)
    {
      assert(sequence.at(idx) != wildcard);
    }
    wildcard_start = std::get<1>(nwr_it) + 1;
  }
  for (size_t idx = wildcard_start; idx < sequence.size(); idx++)
  {
    assert(sequence.at(idx) == wildcard);
  }
}
#endif

static int count_non_wildcard_ranges(const char *progname,
                                     const char *inputfilename)
{
  constexpr const int buf_size = 1 << 14;
  int ret = guess_if_protein_file(progname,inputfilename);
  if (ret == -1)
  {
    return -1;
  }
  if (ret != 0)
  {
    std::cerr << "Usage: " << progname << ": can only handle DNA sequences" 
              << std::endl;
    return -1;
  }
  GttlFpType in_fp = gttl_fp_type_open(inputfilename,"rb");
  if (in_fp == nullptr)
  {
    fprintf(stderr,"%s: Cannot open \"%s\"\n",progname, inputfilename);
    return -1;
  }
  GttlSeqIterator<buf_size> gttl_si(in_fp);
  size_t total_length = 0;
  size_t num_of_sequences = 0;
  size_t non_wildcard_ranges_total_length = 0;
  size_t non_wildcard_ranges_max_length = 0;
  size_t non_wildcard_ranges_number = 0;
  size_t max_sequence_length = 0;
  bool haserr = false;
  try
  {
    for (auto &&si : gttl_si)
    {
      auto sequence = std::get<1>(si);
      total_length += sequence.size();
      max_sequence_length = std::max(max_sequence_length,sequence.size());
      NonWildCardRangeIterator<'N'> nwcr_it(sequence.data(),sequence.size());
      std::vector<std::pair<size_t,size_t>> sequence_ranges 
        = nwcr_it.enumerate();
#ifndef NDEBUG
      verify_non_wildcard_ranges<'N'>(sequence,sequence_ranges);
#endif
      for (auto &&nwr_it : sequence_ranges)
      {
        std::pair<uint64_t,size_t> result;
        const size_t this_length = std::get<1>(nwr_it) -
                                   std::get<0>(nwr_it) + 1;
        non_wildcard_ranges_total_length += this_length;
        non_wildcard_ranges_max_length
          = std::max(non_wildcard_ranges_max_length,this_length);
      }
      non_wildcard_ranges_number += sequence_ranges.size();
      num_of_sequences++;
    }
  }
  catch (std::ios_base::failure const &ext)
  {
    std::cerr << progname << ": " << ext.what() << std::endl;
    haserr = true;
  }
  if (!haserr)
  {
    printf("# num_of_sequences\t%lu\n",num_of_sequences);
    printf("# total_length\t%lu\n",total_length);
    printf("# max_sequence_length\t%lu\n",max_sequence_length);
    printf("# non_wildcard_ranges_number\t%lu\n",non_wildcard_ranges_number);
    printf("# non_wildcard_ranges_total_length\t%lu\n",
              non_wildcard_ranges_total_length);
    printf("# non_wildcard_ranges_max_length\t%lu\n",
              non_wildcard_ranges_max_length);
    printf("# num_of_sequences.bits\t%d\n",gt_required_bits(num_of_sequences));
    printf("# max_sequence_length.bits\t%d\n",
            gt_required_bits<size_t>(max_sequence_length));
    printf("# non_wildcard_ranges_number.bits\t%d\n",
              gt_required_bits<size_t>(non_wildcard_ranges_number));
    printf("# non_wildcard_ranges_max_length.bits\t%d\n",
              gt_required_bits<size_t>(non_wildcard_ranges_max_length));
  }
  gttl_fp_type_close(in_fp);
  return haserr ? -1 : 0;
}

int main(int argc,char *argv[])
{
  const char *progname = argv[0];
  bool haserr = false;
  for (int idx = 1; idx < argc; idx++)
  {
    if (count_non_wildcard_ranges(progname,argv[idx]) != 0)
    {
      haserr = true;
      break;
    }
  }
  return haserr ? EXIT_FAILURE : EXIT_SUCCESS;
}
