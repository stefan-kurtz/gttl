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
#ifndef QGRAMS_HASH2_INVINT_HPP
#define QGRAMS_HASH2_INVINT_HPP
#include "sequences/alphabet.hpp"
#include "sequences/qgrams_rec_hash2_value_iter.hpp"

static int aux_data_invint4(size_t qgram_length)
{
  return static_cast<int>(2 * (qgram_length-1));
}

static uint64_t first_qgram_invint4(const uint8_t *t_qgram,
                                    size_t qgram_length)
{
  uint64_t code = 0, mult= 1;
  for (const uint8_t *t_qgram_ptr = t_qgram + qgram_length - 1;
       t_qgram_ptr >= t_qgram; t_qgram_ptr--)
  {
    code += mult * static_cast<uint64_t>(*t_qgram_ptr);
    mult *= 4;
  }
  return code;
}

static uint64_t first_compl_qgram_invint4(const uint8_t *t_qgram,
                                          size_t qgram_length)
{
  uint64_t code = 0, mult= 1;
  for (const uint8_t *t_qgram_ptr = t_qgram;
       t_qgram_ptr < t_qgram + qgram_length; t_qgram_ptr++)
  {
    code += mult * static_cast<uint64_t>(*t_qgram_ptr);
    mult *= 4;
  }
  return code;
}

static uint64_t next_qgram_invint4(uint8_t old_t_char,
                                   uint64_t integer_code,
                                   uint8_t new_t_char,
                                   int &shift)
{
  integer_code -= (static_cast<uint64_t>(old_t_char) << shift);
  integer_code *= static_cast<uint64_t>(4);
  integer_code += static_cast<uint64_t>(new_t_char);
  return integer_code;
}

static uint64_t next_compl_qgram_invint4(uint8_t old_t_char,
                                         uint64_t integer_code,
                                         uint8_t new_t_char,
                                         int &shift)
{
  static constexpr const uint64_t complement[] = {uint64_t(2),uint64_t(3),
                                                  uint64_t(0),uint64_t(1)};
  integer_code -= complement[old_t_char];
  integer_code /= static_cast<uint64_t>(4);
  integer_code += complement[new_t_char] << shift;
  return integer_code;
}

using InvertibleIntegercode2Iterator4
  = QgramRecHash2ValueIterator<first_qgram_invint4,
                               first_compl_qgram_invint4,
                               int,
                               aux_data_invint4,
                               next_qgram_invint4,
                               next_compl_qgram_invint4>;
#endif