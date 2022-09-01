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
#include <cstdint>
#include <cassert>

class NThashTransformer
{
  // rotate "v" to the left 1 position
  uint64_t rotate_left_1(const uint64_t v) const noexcept
  {
    return (v << 1) | (v >> 63);
  }

  // rotate "v" to the right 1 position
  uint64_t rotate_right_1(const uint64_t v) const noexcept
  {
    return (v >> 1) | (v << 63);
  }

  // swap bit 0 with bit 33 in "v"
  uint64_t swapbits033(const uint64_t v) const noexcept
  {
    uint64_t x = (v ^ (v >> 33)) & uint64_t(1);
    return v ^ (x | (x << 33));
  }

  // swap bit 32 with bit 63 in "v"
  inline uint64_t swapbits3263(const uint64_t v) const noexcept
  {
    uint64_t x = ((v >> 32) ^ (v >> 63)) & 1;
    return v ^ ((x << 32) | (x << 63));
  }

  // 64-bit random seeds corresponding to bases and their complements
  static constexpr const uint64_t nt_hash_seedA = 0x3c8bfbb395c60474;
  static constexpr const uint64_t nt_hash_seedC = 0x3193c18562a02b4c;
  static constexpr const uint64_t nt_hash_seedG = 0x20323ed082572324;
  static constexpr const uint64_t nt_hash_seedT = 0x295549f54be24456;
  static const uint64_t nt_hash_seedN = 0x0000000000000000;

  static constexpr const uint64_t nt_hash_seed_table[] = {
      nt_hash_seedA, nt_hash_seedC, nt_hash_seedG, nt_hash_seedT,
      nt_hash_seedN};

  static constexpr const uint64_t nt_hash_A33r[33] = {
      0x195c60474, 0x12b8c08e9, 0x571811d3,  0xae3023a6,  0x15c60474c,
      0xb8c08e99,  0x171811d32, 0xe3023a65,  0x1c60474ca, 0x18c08e995,
      0x11811d32b, 0x3023a657,  0x60474cae,  0xc08e995c,  0x1811d32b8,
      0x1023a6571, 0x474cae3,   0x8e995c6,   0x11d32b8c,  0x23a65718,
      0x474cae30,  0x8e995c60,  0x11d32b8c0, 0x3a657181,  0x74cae302,
      0xe995c604,  0x1d32b8c08, 0x1a6571811, 0x14cae3023, 0x995c6047,
      0x132b8c08e, 0x6571811d,  0xcae3023a};

  static constexpr const uint64_t nt_hash_A31l[31] = {
      0x3c8bfbb200000000, 0x7917f76400000000, 0xf22feec800000000,
      0xe45fdd9200000000, 0xc8bfbb2600000000, 0x917f764e00000000,
      0x22feec9e00000000, 0x45fdd93c00000000, 0x8bfbb27800000000,
      0x17f764f200000000, 0x2feec9e400000000, 0x5fdd93c800000000,
      0xbfbb279000000000, 0x7f764f2200000000, 0xfeec9e4400000000,
      0xfdd93c8a00000000, 0xfbb2791600000000, 0xf764f22e00000000,
      0xeec9e45e00000000, 0xdd93c8be00000000, 0xbb27917e00000000,
      0x764f22fe00000000, 0xec9e45fc00000000, 0xd93c8bfa00000000,
      0xb27917f600000000, 0x64f22fee00000000, 0xc9e45fdc00000000,
      0x93c8bfba00000000, 0x27917f7600000000, 0x4f22feec00000000,
      0x9e45fdd800000000};

  static constexpr const uint64_t nt_hash_C33r[33] = {
      0x162a02b4c, 0xc5405699,  0x18a80ad32, 0x115015a65, 0x2a02b4cb,
      0x54056996,  0xa80ad32c,  0x15015a658, 0xa02b4cb1,  0x140569962,
      0x80ad32c5,  0x1015a658a, 0x2b4cb15,   0x569962a,   0xad32c54,
      0x15a658a8,  0x2b4cb150,  0x569962a0,  0xad32c540,  0x15a658a80,
      0xb4cb1501,  0x169962a02, 0xd32c5405,  0x1a658a80a, 0x14cb15015,
      0x9962a02b,  0x132c54056, 0x658a80ad,  0xcb15015a,  0x1962a02b4,
      0x12c540569, 0x58a80ad3,  0xb15015a6};

  static constexpr const uint64_t nt_hash_C31l[31] = {
      0x3193c18400000000, 0x6327830800000000, 0xc64f061000000000,
      0x8c9e0c2200000000, 0x193c184600000000, 0x3278308c00000000,
      0x64f0611800000000, 0xc9e0c23000000000, 0x93c1846200000000,
      0x278308c600000000, 0x4f06118c00000000, 0x9e0c231800000000,
      0x3c18463200000000, 0x78308c6400000000, 0xf06118c800000000,
      0xe0c2319200000000, 0xc184632600000000, 0x8308c64e00000000,
      0x6118c9e00000000,  0xc23193c00000000,  0x1846327800000000,
      0x308c64f000000000, 0x6118c9e000000000, 0xc23193c000000000,
      0x8463278200000000, 0x8c64f0600000000,  0x118c9e0c00000000,
      0x23193c1800000000, 0x4632783000000000, 0x8c64f06000000000,
      0x18c9e0c200000000};

  static constexpr const uint64_t nt_hash_G33r[33] = {
      0x82572324,  0x104ae4648, 0x95c8c91,   0x12b91922,  0x25723244,
      0x4ae46488,  0x95c8c910,  0x12b919220, 0x57232441,  0xae464882,
      0x15c8c9104, 0xb9192209,  0x172324412, 0xe4648825,  0x1c8c9104a,
      0x191922095, 0x12324412b, 0x46488257,  0x8c9104ae,  0x11922095c,
      0x324412b9,  0x64882572,  0xc9104ae4,  0x1922095c8, 0x124412b91,
      0x48825723,  0x9104ae46,  0x122095c8c, 0x4412b919,  0x88257232,
      0x1104ae464, 0x2095c8c9,  0x412b9192};

  static constexpr const uint64_t nt_hash_G31l[31] = {
      0x20323ed000000000, 0x40647da000000000, 0x80c8fb4000000000,
      0x191f68200000000,  0x323ed0400000000,  0x647da0800000000,
      0xc8fb41000000000,  0x191f682000000000, 0x323ed04000000000,
      0x647da08000000000, 0xc8fb410000000000, 0x91f6820200000000,
      0x23ed040600000000, 0x47da080c00000000, 0x8fb4101800000000,
      0x1f68203200000000, 0x3ed0406400000000, 0x7da080c800000000,
      0xfb41019000000000, 0xf682032200000000, 0xed04064600000000,
      0xda080c8e00000000, 0xb410191e00000000, 0x6820323e00000000,
      0xd040647c00000000, 0xa080c8fa00000000, 0x410191f600000000,
      0x820323ec00000000, 0x40647da00000000,  0x80c8fb400000000,
      0x10191f6800000000};

  static constexpr const uint64_t nt_hash_T33r[33] = {
      0x14be24456, 0x97c488ad,  0x12f89115a, 0x5f1222b5,  0xbe24456a,
      0x17c488ad4, 0xf89115a9,  0x1f1222b52, 0x1e24456a5, 0x1c488ad4b,
      0x189115a97, 0x11222b52f, 0x24456a5f,  0x488ad4be,  0x9115a97c,
      0x1222b52f8, 0x4456a5f1,  0x88ad4be2,  0x1115a97c4, 0x22b52f89,
      0x456a5f12,  0x8ad4be24,  0x115a97c48, 0x2b52f891,  0x56a5f122,
      0xad4be244,  0x15a97c488, 0xb52f8911,  0x16a5f1222, 0xd4be2445,
      0x1a97c488a, 0x152f89115, 0xa5f1222b};

  static constexpr const uint64_t nt_hash_T31l[31] = {
      0x295549f400000000, 0x52aa93e800000000, 0xa55527d000000000,
      0x4aaa4fa200000000, 0x95549f4400000000, 0x2aa93e8a00000000,
      0x55527d1400000000, 0xaaa4fa2800000000, 0x5549f45200000000,
      0xaa93e8a400000000, 0x5527d14a00000000, 0xaa4fa29400000000,
      0x549f452a00000000, 0xa93e8a5400000000, 0x527d14aa00000000,
      0xa4fa295400000000, 0x49f452aa00000000, 0x93e8a55400000000,
      0x27d14aaa00000000, 0x4fa2955400000000, 0x9f452aa800000000,
      0x3e8a555200000000, 0x7d14aaa400000000, 0xfa29554800000000,
      0xf452aa9200000000, 0xe8a5552600000000, 0xd14aaa4e00000000,
      0xa295549e00000000, 0x452aa93e00000000, 0x8a55527c00000000,
      0x14aaa4fa00000000};

  static constexpr const uint64_t nt_hash_N33r[33] = {
      nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN,
      nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN,
      nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN,
      nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN,
      nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN,
      nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN,
      nt_hash_seedN, nt_hash_seedN, nt_hash_seedN};

  static constexpr const uint64_t nt_hash_N31l[31] = {
      nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN,
      nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN,
      nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN,
      nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN,
      nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN,
      nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN, nt_hash_seedN,
      nt_hash_seedN};

  static constexpr const uint64_t *msTab31l[] = {
      nt_hash_A31l, nt_hash_C31l, nt_hash_G31l, nt_hash_T31l, nt_hash_N31l};
  static constexpr const uint64_t *msTab33r[] = {
      nt_hash_A33r, nt_hash_C33r, nt_hash_G33r, nt_hash_T33r, nt_hash_N33r};
  uint64_t msTab31l_33r_or[4 + 1];

 public:
  NThashTransformer(size_t qgram_length)
  {
    for (size_t idx = 0; idx < sizeof msTab31l_33r_or/sizeof msTab31l_33r_or[0];
         idx++)
    {
      msTab31l_33r_or[idx]
        = (msTab31l[idx][qgram_length % 31] | msTab33r[idx][qgram_length % 33]);
    }
  }
  uint64_t first_hash_value_get(const uint8_t *kmerSeq, size_t qgram_length)
    const noexcept
  {
    uint64_t hVal = 0;
    for (size_t idx = 0; idx < qgram_length; idx++)
    {
      hVal = rotate_left_1(hVal);
      hVal = swapbits033(hVal);
      assert(kmerSeq[idx] <= 4);
      hVal ^= nt_hash_seed_table[kmerSeq[idx]];
    }
    return hVal;
  }
  // forward-strand ntHash for sliding k-mers
  uint64_t next_hash_value_get(uint8_t charOut, uint64_t fhVal,
                               uint8_t charIn) const noexcept
  {
    assert(charIn <= 4 && charOut <= 4);
    uint64_t hVal = rotate_left_1(fhVal);
    hVal = swapbits033(hVal);
    hVal ^= nt_hash_seed_table[charIn];
    hVal ^= msTab31l_33r_or[charOut];
    return hVal;
  }
  // reverse-strand ntHash for sliding k-mers
  uint64_t next_compl_hash_value_get(uint8_t compl_charOut,
                                     uint64_t rhVal,
                                     uint8_t compl_charIn)
                                     const noexcept
  {
    assert(compl_charIn < 4 && compl_charOut < 4);
    uint64_t hVal = rhVal ^ msTab31l_33r_or[compl_charIn];
    hVal ^= nt_hash_seed_table[compl_charOut];
    hVal = rotate_right_1(hVal);
    hVal = swapbits3263(hVal);
    return hVal;
  }
};
#endif
