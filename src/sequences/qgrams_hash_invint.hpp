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
#ifndef QGRAMS_HASH_INVINT_HPP
#define QGRAMS_HASH_INVINT_HPP
#include "sequences/alphabet.hpp"
#include "sequences/qgrams_hash_rec_iter.hpp"

static uint64_t next_qgram_integer_code_asize4(uint8_t old_t_char,
                                               uint64_t integer_code,
                                               uint8_t new_t_char,
                                               int &shift)
{
  integer_code -= (static_cast<uint64_t>(old_t_char) << shift);
  integer_code *= static_cast<uint64_t>(4);
  integer_code += static_cast<uint64_t>(new_t_char);
  return integer_code;
}

static int aux_data_integer_code_asize4(size_t qgram_length)
{
  return static_cast<int>(2 * (qgram_length-1));
}

static uint64_t next_qgram_integer_code_asize20(uint8_t old_t_char,
                                                uint64_t integer_code,
                                                uint8_t new_t_char,
                                                uint64_t &msb_weight)
{
  integer_code -= static_cast<uint64_t>(old_t_char) * msb_weight;
  integer_code *= static_cast<uint64_t>(20);
  integer_code += static_cast<uint64_t>(new_t_char);
  return integer_code;
}

template<size_t alpha_size>
static uint64_t first_qgram_integer_code(const uint8_t *t_qgram,
                                         size_t qgram_length)
{
  uint64_t code = 0, mult= 1;
  for (const uint8_t *t_qgram_ptr = t_qgram + qgram_length - 1;
       t_qgram_ptr >= t_qgram; t_qgram_ptr--)
  {
    code += mult * static_cast<uint64_t>(*t_qgram_ptr);
    mult *= alpha_size;
  }
  return code;
}

static uint64_t aux_data_integer_code_asize20(size_t qgram_length)
{
  return static_cast<uint64_t>(std::pow(static_cast<size_t>(20),
                                        qgram_length-1));
}

using InvertibleIntegercodeIterator4
  = QgramRecHashValueIterator<alphabet::nucleotides_upper_lower,
                              4,
                              first_qgram_integer_code<4>,
                              int,
                              aux_data_integer_code_asize4,
                              next_qgram_integer_code_asize4>;
using InvertibleIntegercodeIterator4_Wildcard2_a
  = QgramRecHashValueIterator<alphabet::nucleotides_upper_lower,
                              0,
                              first_qgram_integer_code<4>,
                              int,
                              aux_data_integer_code_asize4,
                              next_qgram_integer_code_asize4>;
using InvertibleIntegercodeIterator20
  = QgramRecHashValueIterator<alphabet::amino_acids,
                              20,
                              first_qgram_integer_code<20>,
                              uint64_t,
                              aux_data_integer_code_asize20,
                              next_qgram_integer_code_asize20>;
#endif
