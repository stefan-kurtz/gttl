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
#ifndef GUESS_IF_PROTEIN_SEQ_HPP
#define GUESS_IF_PROTEIN_SEQ_HPP
#include <cstddef>
#include <cstring>
#include <climits>
#include <algorithm>
#include <ios>
#include <vector>
#include <string>
#include "sequences/gttl_fasta_generator.hpp"
#include "utilities/gttl_file_open.hpp"

/* checks first chars of sequence, returns true if sequence is
   definitely a protein sequence */

inline constexpr size_t GUESS_SIZE_TO_DECIDE = 1000;

inline bool guess_if_protein_sequence(const char *sequence,size_t seqlen)
{
  const char *const protein_only_characters = "LIFEQPXZ";
  bool for_protein_only_lookup[UCHAR_MAX+1] = {false};
  for (size_t idx = 0; idx < strlen(protein_only_characters); idx++)
  {
    for_protein_only_lookup[static_cast<int>(protein_only_characters[idx])] =
                            true;
  }
  for (size_t idx = 0; idx < std::min(GUESS_SIZE_TO_DECIDE,seqlen);
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

template<class ReaderClass>
bool guess_if_protein_file_generic(const char *filename)
{
  const GttlFpType in_fp = gttl_fp_type_open(filename, "rb");

  if (in_fp == nullptr)
  {
    throw std::ios_base::failure(": cannot open file");
  }
  ReaderClass reader(in_fp);
  size_t total_length = 0;
  bool decided_if_protein = false;
  for (auto &&si : reader)
  {
    auto sequence = si->sequence_get();
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
  return decided_if_protein;
}

inline bool guess_if_protein_file(const char *filename)
{
  constexpr const int buf_size = 1 << 14;
  return guess_if_protein_file_generic<GttlFastAGenerator<buf_size>>(filename);
}

inline bool guess_if_protein_file(const std::vector<std::string> &inputfiles)
{
  for (const std::string& file : inputfiles)
  {
    if (guess_if_protein_file(file.c_str())) return true;
  }
  return false;
}

template<class MultiseqClass>
bool guess_if_protein_multiseq(const MultiseqClass *multiseq)
{
  const size_t sequences_number = multiseq->sequences_number_get();
  bool is_protein_sequence = false;
  size_t sequences_total_length = 0;
  for (size_t seqnum = 0; seqnum < sequences_number; seqnum++)
  {
    const char *const seqptr = multiseq->sequence_ptr_get(seqnum);
    const size_t seqlen = multiseq->sequence_length_get(seqnum);
    sequences_total_length += seqlen;
    if (guess_if_protein_sequence(seqptr,seqlen))
    {
      is_protein_sequence = true;
      break;
    }
    if (sequences_total_length >= GUESS_SIZE_TO_DECIDE)
    {
      break;
    }
  }
  return is_protein_sequence;
}
#endif
