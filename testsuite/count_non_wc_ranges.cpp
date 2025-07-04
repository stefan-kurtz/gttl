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
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <vector>
#include <utility>
#include "sequences/gttl_fasta_generator.hpp"
#include "utilities/gttl_file_open.hpp"
#include "utilities/mathsupport.hpp"
#include "sequences/guess_if_protein_seq.hpp"
#include "sequences/non_wildcard_ranges.hpp"

#ifndef NDEBUG
template<char wildcard>
static void verify_non_wildcard_ranges(const std::string_view &sequence,
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

static void count_non_wildcard_ranges(const char *inputfilename)
{
  constexpr const int buf_size = 1 << 14;
  const bool is_protein = guess_if_protein_file(inputfilename);

  if (is_protein)
  {
    throw std::runtime_error(": can only handle DNA sequences");
    /* check_err.py checked */
  }
  const GttlFpType in_fp = gttl_fp_type_open(inputfilename, "rb");
  if (in_fp == nullptr)
  {
    throw std::runtime_error(": cannot open file");
    /* check_err.py checked */
  }
  GttlFastAGenerator<buf_size> gttl_si(in_fp);
  size_t total_length = 0;
  size_t num_of_sequences = 0;
  size_t non_wildcard_ranges_total_length = 0;
  size_t non_wildcard_ranges_max_length = 0;
  size_t non_wildcard_ranges_number = 0;
  size_t max_sequence_length = 0;
  for (auto &&si : gttl_si)
  {
    auto sequence = si->sequence_get();
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
      const size_t this_length = std::get<1>(nwr_it) -
                                 std::get<0>(nwr_it) + 1;
      non_wildcard_ranges_total_length += this_length;
      non_wildcard_ranges_max_length
        = std::max(non_wildcard_ranges_max_length,this_length);
    }
    non_wildcard_ranges_number += sequence_ranges.size();
    num_of_sequences++;
  }
  printf("# filename\t%s\n",inputfilename);
  printf("# num_of_sequences\t%zu\n",num_of_sequences);
  printf("# total_length\t%zu\n",total_length);
  printf("# max_sequence_length\t%zu\n",max_sequence_length);
  printf("# non_wildcard_ranges_number\t%zu\n",non_wildcard_ranges_number);
  printf("# non_wildcard_ranges_total_length\t%zu\n",
            non_wildcard_ranges_total_length);
  printf("# non_wildcard_ranges_max_length\t%zu\n",
            non_wildcard_ranges_max_length);
  printf("# num_of_sequences.bits\t%d\n",gttl_required_bits(num_of_sequences));
  printf("# max_sequence_length.bits\t%d\n",
          gttl_required_bits<size_t>(max_sequence_length));
  printf("# non_wildcard_ranges_number.bits\t%d\n",
            gttl_required_bits<size_t>(non_wildcard_ranges_number));
  printf("# non_wildcard_ranges_max_length.bits\t%d\n",
            gttl_required_bits<size_t>(non_wildcard_ranges_max_length));
}

int main(int argc,char *argv[])
{
  const char *const progname = argv[0];
  bool haserr = false;
  for (int idx = 1; idx < argc; idx++)
  {
    try
    {
      count_non_wildcard_ranges(argv[idx]);
    }
    catch (const std::exception &msg)
    {
      std::cerr << progname << ": file \"" << argv[idx] << "\"" << msg.what()
                << '\n';
      haserr = true;
      break;
    }
  }
  return haserr ? EXIT_FAILURE : EXIT_SUCCESS;
}
