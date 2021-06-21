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
#ifndef GUESS_IF_PROTEIN_SEQ_HPP_
#define GUESS_IF_PROTEIN_SEQ_HPP_
#include <cstdbool>
#include <cstddef>
#include <cstring>
#include <climits>
#include <algorithm>
#include <iostream>
#include "sequences/gttl_seq_iterator.hpp"

/* checks first chars of sequence, returns true if sequence is
   definitely a protein sequence */

#define GUESS_SIZE_TO_DECIDE 1000

bool guess_if_protein_sequence(const char *sequence,size_t seqlen)
{
  const char *protein_only_characters = "LIFEQPXZ";
  bool for_protein_only_lookup[UCHAR_MAX+1] = {false};
  for (size_t idx = 0; idx < strlen(protein_only_characters); idx++)
  {
    for_protein_only_lookup[static_cast<int>(protein_only_characters[idx])] =
                            true;
  }
  for (size_t idx = 0; idx < std::min((size_t) GUESS_SIZE_TO_DECIDE,seqlen);
       idx++)
  {
    if (for_protein_only_lookup[static_cast<int>(sequence[idx])])
    {
      // no need to check more bases, so break for loop
      return true;
    }
  }
  return false;
}

int guess_if_protein_file(const char *progname,const char *filename)
{
  const int buf_size = 1 << 14;
  GttlFpType in_fp = gttl_fp_type_open(filename,"rb");

  if (in_fp == nullptr)
  {
    fprintf(stderr,"%s: Cannot open \"%s\"\n",progname, filename);
    return -1;
  }
  GttlSeqIterator<buf_size> gttl_si(in_fp);
  size_t total_length = 0;
  bool haserr = false;
  bool decided_if_protein = false;
  try
  {
    for (auto &&si : gttl_si)
    {
      auto sequence = std::get<1>(si);
      if (guess_if_protein_sequence(sequence.data(),sequence.size()))
      {
        decided_if_protein = true;
        break;
      }
      total_length += sequence.size();
      if (total_length >= GUESS_SIZE_TO_DECIDE)
      {
        break;
      }
    }
  }
  catch (std::ios_base::failure const &ext)
  {
    std::cerr << progname << ": " << ext.what() << std::endl;
    haserr = true;
  }
  gttl_fp_type_close(in_fp);
  if (haserr)
  {
    return -1;
  }
  if (decided_if_protein) /* have found evidence of protein in
                             GUESS_SIZE_TO_DECIDE first characters */
  {
    return 1; /* it is a protein */
  }
  return 0;
}
#endif
