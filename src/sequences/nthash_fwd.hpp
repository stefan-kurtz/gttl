/*
MIT License

Copyright (c) 2018 Hamid Mohamadi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
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
#ifndef NTHASH_FWD_HPP
#define NTHASH_FWD_HPP
#include <array>
#include <cstddef>
#include <cstdint>
#include <cassert>
#include <utility>
#include "sequences/complement_uint8.hpp"
#include "sequences/nthash_rotation_tables.hpp"

class NThashTransformer
{
  // rotate "v" to the left 1 position
  [[nodiscard]] static uint64_t rotate_left_1(const uint64_t v) noexcept
  {
    return (v << 1) | (v >> 63);
  }

  // rotate "v" to the right 1 position
  [[nodiscard]] static uint64_t rotate_right_1(const uint64_t v) noexcept
  {
    return (v >> 1) | (v << 63);
  }

  // swap bit 0 with bit 33 in "v"
  [[nodiscard]] static uint64_t swapbits033(const uint64_t v) noexcept
  {
    const uint64_t x = (v ^ (v >> 33)) & uint64_t{1};
    return v ^ (x | (x << 33));
  }

  // swap bit 32 with bit 63 in "v"
  [[nodiscard]] static uint64_t swapbits3263(const uint64_t v) noexcept
  {
    const uint64_t x = ((v >> 32) ^ (v >> 63)) & 1;
    return v ^ ((x << 32) | (x << 63));
  }

  // 64-bit random seeds corresponding to bases and their complements
  static constexpr const uint64_t nt_hash_seedA = 0x3c8bfbb395c60474;
  static constexpr const uint64_t nt_hash_seedC = 0x3193c18562a02b4c;
  static constexpr const uint64_t nt_hash_seedG = 0x20323ed082572324;
  static constexpr const uint64_t nt_hash_seedT = 0x295549f54be24456;

  static constexpr const uint64_t nt_hash_seed_table[] = {
      nt_hash_seedA, nt_hash_seedC, nt_hash_seedG, nt_hash_seedT};

  static constexpr const auto nt_hash_A33r = generate_33r_table(nt_hash_seedA);
  static constexpr const auto nt_hash_C33r = generate_33r_table(nt_hash_seedC);
  static constexpr const auto nt_hash_G33r = generate_33r_table(nt_hash_seedG);
  static constexpr const auto nt_hash_T33r = generate_33r_table(nt_hash_seedT);
  static constexpr const auto nt_hash_A31l = generate_31l_table(nt_hash_seedA);
  static constexpr const auto nt_hash_C31l = generate_31l_table(nt_hash_seedC);
  static constexpr const auto nt_hash_G31l = generate_31l_table(nt_hash_seedG);
  static constexpr const auto nt_hash_T31l = generate_31l_table(nt_hash_seedT);

  static constexpr const std::array<std::array<uint64_t, 31>, 4> msTab31l = {
    nt_hash_A31l, nt_hash_C31l, nt_hash_G31l, nt_hash_T31l};
  static constexpr const std::array<std::array<uint64_t, 33>, 4> msTab33r = {
    nt_hash_A33r, nt_hash_C33r, nt_hash_G33r, nt_hash_T33r};

  uint64_t msTab31l_33r_or[4] = {};

 public:
  static constexpr const bool possible_false_positive_matches = true;
  explicit NThashTransformer(size_t qgram_length)
  {
    for (size_t idx = 0; idx < sizeof msTab31l_33r_or/sizeof msTab31l_33r_or[0];
         idx++)
    {
      msTab31l_33r_or[idx]
        = (msTab31l[idx][qgram_length % 31] | msTab33r[idx][qgram_length % 33]);
    }
  }
  static uint64_t first_fwd_hash_value_get(const uint8_t *t_sequence,
                                    size_t qgram_length)  noexcept
  {
    uint64_t hVal = 0;
    for (size_t idx = 0; idx < qgram_length; idx++)
    {
      hVal = rotate_left_1(hVal);
      hVal = swapbits033(hVal);
      assert(t_sequence[idx] < uint8_t(4));
      hVal ^= nt_hash_seed_table[t_sequence[idx]];
    }
    return hVal;
  }
  static std::pair<uint64_t,uint64_t> first_hash_value_pair_get(
                                    const uint8_t *t_qgram,
                                    size_t qgram_length)  noexcept
  {
    const uint64_t fwd_code = first_fwd_hash_value_get(t_qgram,qgram_length);
    uint64_t rev_compl_code = 0;
    for (size_t idx = qgram_length; idx > 0; /*Nothing*/)
    {
      idx--;
      rev_compl_code = rotate_left_1(rev_compl_code);
      rev_compl_code = swapbits033(rev_compl_code);
      const uint8_t complement_cc = complement_uint8(t_qgram[idx]);
      rev_compl_code ^= nt_hash_seed_table[static_cast<size_t>(complement_cc)];
    }
    return std::make_pair(fwd_code,rev_compl_code);
  }

  // forward-strand ntHash for sliding k-mers
  [[nodiscard]] uint64_t next_hash_value_get(uint8_t charOut, uint64_t fhVal,
                                             uint8_t charIn) const noexcept
  {
    assert(charIn < uint8_t(4) && charOut < uint8_t(4));
    uint64_t hVal = rotate_left_1(fhVal);
    hVal = swapbits033(hVal);
    hVal ^= nt_hash_seed_table[charIn];
    hVal ^= msTab31l_33r_or[charOut];
    return hVal;
  }
  // reverse-strand ntHash for sliding k-mers
  [[nodiscard]] uint64_t next_compl_hash_value_get(uint8_t compl_charOut,
                                     uint64_t rhVal,
                                     uint8_t compl_charIn)
                                     const noexcept
  {
    assert(compl_charIn < uint8_t(4) && compl_charOut < uint8_t(4));
    uint64_t hVal = rhVal ^ msTab31l_33r_or[compl_charIn];
    hVal ^= nt_hash_seed_table[compl_charOut];
    hVal = rotate_right_1(hVal);
    hVal = swapbits3263(hVal);
    return hVal;
  }
};
#endif
