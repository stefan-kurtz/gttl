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
#include "utilities/mathsupport.hpp"
#include "sequences/gttl_seq_iterator.hpp"
#include "sequences/guess_if_protein_seq.hpp"
#include "sequences/char_range.hpp"

static void display_char_ranges(const char *inputfilename)
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
  size_t non_wildcard_ranges_total_length = 0;
  try /* need this, as the catch needs to close the file pointer
         to prevent a memory leak */
  {
    size_t seqnum = 0;
    for (auto &&si : gttl_si)
    {
      auto sequence = std::get<1>(si);
      static constexpr const char nucleotides[] = "ACGTacgt";
      GttlCharRange<nucleotides> ranger(sequence.data(),sequence.size());
      for (auto &&it : ranger)
      {
        std::cout << seqnum << "\t" << std::get<0>(it)
                  << "\t" << std::get<1>(it) << std::endl;
        non_wildcard_ranges_total_length += std::get<1>(it);
      }
      seqnum++;
    }
  }
  catch (std::string &msg)
  {
    gttl_fp_type_close(in_fp);
    throw msg;
  }
  gttl_fp_type_close(in_fp);
  std::cout << "# non_wildcard_ranges_total_length\t" <<
                non_wildcard_ranges_total_length << std::endl;
}

int main(int argc,char *argv[])
{
  const char *progname = argv[0];
  bool haserr = false;
  for (int idx = 1; idx < argc; idx++)
  {
    try
    {
      display_char_ranges(argv[idx]);
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
