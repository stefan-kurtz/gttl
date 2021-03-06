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
#include "sequences/qgrams_rec_hash_value_iter.hpp"
#include "sequences/qgrams_rec_hash2_value_iter.hpp"

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

class InvertibleIntegercodeTransformer4
{
  static constexpr const uint64_t invint4_complement_table[]
    = {uint64_t(2),uint64_t(3), uint64_t(0),uint64_t(1)};
  int shift;
  public:
  static constexpr const char nucleotides_upper[] = "ACTG";
  InvertibleIntegercodeTransformer4(size_t qgram_length)
    : shift(static_cast<int>(2 * (qgram_length-1)))
  {}
  uint64_t first_hash_value_get(const uint8_t *sequence, size_t qgram_length)
    const noexcept
  {
    return first_qgram_integer_code<4>(sequence,qgram_length);
  }
  uint64_t next_hash_value_get(uint8_t old_t_char,
                               uint64_t integer_code,
                               uint8_t new_t_char) const noexcept
  {
    integer_code -= (static_cast<uint64_t>(old_t_char) << shift);
    integer_code *= static_cast<uint64_t>(4);
    integer_code += static_cast<uint64_t>(new_t_char);
    return integer_code;
  }
  uint64_t first_compl_hash_value_get(const uint8_t *t_sequence,
                                      size_t qgram_length) const noexcept
  {
    uint64_t code = 0, mult= 1;
    for (const uint8_t *t_qgram_ptr = t_sequence;
         t_qgram_ptr < t_sequence + qgram_length; t_qgram_ptr++)
    {
      code += mult * invint4_complement_table[*t_qgram_ptr];
      mult *= 4;
    }
    return code;
  }
  uint64_t next_compl_hash_value_get(uint8_t old_t_char,
                                     uint64_t integer_code,
                                     uint8_t new_t_char) const noexcept
  {
    integer_code -= invint4_complement_table[old_t_char];
    integer_code /= static_cast<uint64_t>(4);
    integer_code += (invint4_complement_table[new_t_char] << shift);
    return integer_code;
  }
};

using InvertibleIntegercodeIterator4
  = QgramRecHashValueIterator<alphabet::nucleotides_upper_lower,
                              4,
                              InvertibleIntegercodeTransformer4>;

using InvertibleIntegercodeIterator4_Wildcard2_a
  = QgramRecHashValueIterator<alphabet::nucleotides_upper_lower,
                              0,
                              InvertibleIntegercodeTransformer4>;

using InvertibleIntegercode2Iterator4
  = QgramRecHash2ValueIterator<InvertibleIntegercodeTransformer4>;

class InvertibleIntegercodeTransformer20
{
  uint64_t msb_weight;
  public:
  InvertibleIntegercodeTransformer20(size_t qgram_length)
    : msb_weight(static_cast<uint64_t>(std::pow(static_cast<size_t>(20),
                                       qgram_length-1)))
  {}
  uint64_t first_hash_value_get(const uint8_t *sequence, size_t qgram_length)
    const noexcept
  {
    return first_qgram_integer_code<20>(sequence,qgram_length);
  }
  uint64_t next_hash_value_get(uint8_t old_t_char,
                               uint64_t integer_code,
                               uint8_t new_t_char) const noexcept
  {
    integer_code -= static_cast<uint64_t>(old_t_char) * msb_weight;
    integer_code *= static_cast<uint64_t>(20);
    integer_code += static_cast<uint64_t>(new_t_char);
    return integer_code;
  }
};

using InvertibleIntegercodeIterator20
  = QgramRecHashValueIterator<alphabet::amino_acids,
                              20,
                              InvertibleIntegercodeTransformer20>;

#endif
