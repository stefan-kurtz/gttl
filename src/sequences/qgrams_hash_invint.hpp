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
#include <cstddef>
#include <cstdint>
#include <utility>
#include <cassert>
#include <cmath>
#include "sequences/alphabet.hpp"
#include "sequences/qgrams_rec_hash_value_fwd_iter.hpp"
#include "sequences/qgrams_rec_hash_value_iter.hpp"

template<size_t alpha_size>
static uint64_t first_fwd_qgram_integer_code(const uint8_t *t_qgram,
                                             size_t qgram_length)
{
  uint64_t code = 0;
  uint64_t mult= 1;
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
  int shift;
  public:
  static constexpr const bool possible_false_positive_matches = false;
  InvertibleIntegercodeTransformer4(size_t qgram_length)
    : shift(static_cast<int>(2 * (qgram_length-1)))
  {}
  uint64_t first_fwd_hash_value_get(const uint8_t *sequence,
                                    size_t qgram_length) const noexcept
  {
    return first_fwd_qgram_integer_code<4>(sequence,qgram_length);
  }
  std::pair<uint64_t,uint64_t> first_hash_value_pair_get(
                                    const uint8_t *t_qgram,
                                    size_t qgram_length) const noexcept
  {
    uint64_t fwd_code = 0;
    uint64_t fwd_mult = 1;
    uint64_t rev_compl_code = 0;
    uint64_t rev_compl_mult= uint64_t(1) << (2 * (qgram_length - 1));
    for (const uint8_t *t_qgram_ptr = t_qgram + qgram_length - 1;
         t_qgram_ptr >= t_qgram; t_qgram_ptr--)
    {
      const uint64_t cc = static_cast<uint64_t>(*t_qgram_ptr);
      fwd_code += fwd_mult * cc;
      fwd_mult *= uint64_t(4);
      assert(cc < uint64_t(4));
      rev_compl_code += rev_compl_mult * (uint64_t(3) - cc);
      rev_compl_mult /= uint64_t(4);
    }
    return std::make_pair(fwd_code,rev_compl_code);
  }
  [[nodiscard]] uint64_t
  next_hash_value_get(uint8_t old_t_char,
                      uint64_t integer_code,
                      uint8_t new_t_char) const noexcept
  {
    integer_code -= (static_cast<uint64_t>(old_t_char) << shift);
    integer_code *= uint64_t(4);
    integer_code += static_cast<uint64_t>(new_t_char);
    return integer_code;
  }
  [[nodiscard]] uint64_t
  next_compl_hash_value_get(uint8_t compl_old_t_char,
                            uint64_t integer_code,
                            uint8_t compl_new_t_char) const noexcept
  {
    assert(integer_code >= static_cast<uint64_t>(compl_old_t_char));
    integer_code -= compl_old_t_char;
    assert(integer_code % uint64_t(4) == 0);
    integer_code /= uint64_t(4);
    integer_code += (static_cast<uint64_t>(compl_new_t_char) << shift);
    return integer_code;
  }
};

using InvertibleIntegercodeIterator4
  = QgramRecHashValueFwdIterator<alphabet::nucleotides_upper_lower,
                                 4,
                                 InvertibleIntegercodeTransformer4,
                                 char>;

using InvertibleIntegercodeIterator4_Wildcard2_a
  = QgramRecHashValueFwdIterator<alphabet::nucleotides_upper_lower,
                                 0,
                                 InvertibleIntegercodeTransformer4,
                                 char>;

using InvertibleIntegercode2Iterator4
  = QgramRecHashValueIterator<alphabet::nucleotides_upper_lower,
                              4,
                              InvertibleIntegercodeTransformer4>;

class InvertibleIntegercodeTransformer20
{
  uint64_t msb_weight;
  public:
  static constexpr const bool possible_false_positive_matches = false;
  InvertibleIntegercodeTransformer20(size_t qgram_length)
    : msb_weight(static_cast<uint64_t>(std::pow(static_cast<size_t>(20),
                                       qgram_length-1)))
  {}
  uint64_t first_fwd_hash_value_get(const uint8_t *sequence,
                                    size_t qgram_length)
    const noexcept
  {
    return first_fwd_qgram_integer_code<20>(sequence,qgram_length);
  }
  [[nodiscard]] uint64_t
  next_hash_value_get(uint8_t old_t_char,
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
  = QgramRecHashValueFwdIterator<alphabet::amino_acids,
                                 20,
                                 InvertibleIntegercodeTransformer20,
                                 char>;
#endif
