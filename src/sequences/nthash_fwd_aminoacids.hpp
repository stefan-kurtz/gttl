#ifndef NTHASH_FWD_AMINOACIDS_HPP
#define NTHASH_FWD_AMINOACIDS_HPP

#include <cassert>
#include <cstddef>
#include <cstdint>
#include "sequences/complement_uint8.hpp"
#include "sequences/nthash_rotation_tables.hpp"
class NtHashAminoacidsTransformer
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
    uint64_t x = (v ^ (v >> 33)) & uint64_t{1};
    return v ^ (x | (x << 33));
  }

  // swap bit 32 with bit 63 in "v"
  [[nodiscard]] static uint64_t swapbits3263(const uint64_t v) noexcept
  {
    uint64_t x = ((v >> 32) ^ (v >> 63)) & 1;
    return v ^ ((x << 32) | (x << 63));
  }

  static constexpr const uint64_t nt_hash_seedA = 0xE6CC5634828DE226;
  static constexpr const uint64_t nt_hash_seedC = 0xC3916818AB69EAEF;
  static constexpr const uint64_t nt_hash_seedD = 0x708B967AB82D90DA;
  static constexpr const uint64_t nt_hash_seedE = 0x97EE132403084238;
  static constexpr const uint64_t nt_hash_seedF = 0x54B5F4F96AA2B4C6;
  static constexpr const uint64_t nt_hash_seedG = 0xA16E4815AEB82F9C;
  static constexpr const uint64_t nt_hash_seedH = 0xD45994107632BA9E;
  static constexpr const uint64_t nt_hash_seedI = 0x52CA8A666024FEB4;
  static constexpr const uint64_t nt_hash_seedK = 0x4B7BE4E0445E0FB4;
  static constexpr const uint64_t nt_hash_seedL = 0x10936D371508CFBB;
  static constexpr const uint64_t nt_hash_seedM = 0x9845685EF8E90AC2;
  static constexpr const uint64_t nt_hash_seedN = 0xBC57ED4598240F30;
  static constexpr const uint64_t nt_hash_seedP = 0x4ADD75731C811286;
  static constexpr const uint64_t nt_hash_seedQ = 0xF4EADE76635CE527;
  static constexpr const uint64_t nt_hash_seedR = 0x725F7D219B2873BC;
  static constexpr const uint64_t nt_hash_seedS = 0x2A90195B4166C058;
  static constexpr const uint64_t nt_hash_seedT = 0xEABF48D47DA98AED;
  static constexpr const uint64_t nt_hash_seedV = 0x7BE88554E92B4BF9;
  static constexpr const uint64_t nt_hash_seedW = 0x211DFFCF403D167C;
  static constexpr const uint64_t nt_hash_seedY = 0x8454A6A3B239AE9E;

  static constexpr const uint64_t nt_hash_seed_table[] = {
    nt_hash_seedA, nt_hash_seedC, nt_hash_seedD, nt_hash_seedE, nt_hash_seedF,
    nt_hash_seedG, nt_hash_seedH, nt_hash_seedI, nt_hash_seedK, nt_hash_seedL,
    nt_hash_seedM, nt_hash_seedN, nt_hash_seedP, nt_hash_seedQ, nt_hash_seedR,
    nt_hash_seedS, nt_hash_seedT, nt_hash_seedV, nt_hash_seedW, nt_hash_seedY
  };

  static constexpr const auto nt_hash_A33r = generate_33r_table(nt_hash_seedA);
  static constexpr const auto nt_hash_C33r = generate_33r_table(nt_hash_seedC);
  static constexpr const auto nt_hash_D33r = generate_33r_table(nt_hash_seedD);
  static constexpr const auto nt_hash_E33r = generate_33r_table(nt_hash_seedE);
  static constexpr const auto nt_hash_F33r = generate_33r_table(nt_hash_seedF);
  static constexpr const auto nt_hash_G33r = generate_33r_table(nt_hash_seedG);
  static constexpr const auto nt_hash_H33r = generate_33r_table(nt_hash_seedH);
  static constexpr const auto nt_hash_I33r = generate_33r_table(nt_hash_seedI);
  static constexpr const auto nt_hash_K33r = generate_33r_table(nt_hash_seedK);
  static constexpr const auto nt_hash_L33r = generate_33r_table(nt_hash_seedL);
  static constexpr const auto nt_hash_M33r = generate_33r_table(nt_hash_seedM);
  static constexpr const auto nt_hash_N33r = generate_33r_table(nt_hash_seedN);
  static constexpr const auto nt_hash_P33r = generate_33r_table(nt_hash_seedP);
  static constexpr const auto nt_hash_Q33r = generate_33r_table(nt_hash_seedQ);
  static constexpr const auto nt_hash_R33r = generate_33r_table(nt_hash_seedR);
  static constexpr const auto nt_hash_S33r = generate_33r_table(nt_hash_seedS);
  static constexpr const auto nt_hash_T33r = generate_33r_table(nt_hash_seedT);
  static constexpr const auto nt_hash_V33r = generate_33r_table(nt_hash_seedV);
  static constexpr const auto nt_hash_W33r = generate_33r_table(nt_hash_seedW);
  static constexpr const auto nt_hash_Y33r = generate_33r_table(nt_hash_seedY);

  static constexpr const auto nt_hash_A31l = generate_31l_table(nt_hash_seedA);
  static constexpr const auto nt_hash_C31l = generate_31l_table(nt_hash_seedC);
  static constexpr const auto nt_hash_D31l = generate_31l_table(nt_hash_seedD);
  static constexpr const auto nt_hash_E31l = generate_31l_table(nt_hash_seedE);
  static constexpr const auto nt_hash_F31l = generate_31l_table(nt_hash_seedF);
  static constexpr const auto nt_hash_G31l = generate_31l_table(nt_hash_seedG);
  static constexpr const auto nt_hash_H31l = generate_31l_table(nt_hash_seedH);
  static constexpr const auto nt_hash_I31l = generate_31l_table(nt_hash_seedI);
  static constexpr const auto nt_hash_K31l = generate_31l_table(nt_hash_seedK);
  static constexpr const auto nt_hash_L31l = generate_31l_table(nt_hash_seedL);
  static constexpr const auto nt_hash_M31l = generate_31l_table(nt_hash_seedM);
  static constexpr const auto nt_hash_N31l = generate_31l_table(nt_hash_seedN);
  static constexpr const auto nt_hash_P31l = generate_31l_table(nt_hash_seedP);
  static constexpr const auto nt_hash_Q31l = generate_31l_table(nt_hash_seedQ);
  static constexpr const auto nt_hash_R31l = generate_31l_table(nt_hash_seedR);
  static constexpr const auto nt_hash_S31l = generate_31l_table(nt_hash_seedS);
  static constexpr const auto nt_hash_T31l = generate_31l_table(nt_hash_seedT);
  static constexpr const auto nt_hash_V31l = generate_31l_table(nt_hash_seedV);
  static constexpr const auto nt_hash_W31l = generate_31l_table(nt_hash_seedW);
  static constexpr const auto nt_hash_Y31l = generate_31l_table(nt_hash_seedY);

  static constexpr const std::array<std::array<uint64_t, 31>, 20> msTab31l = {
    nt_hash_A31l, nt_hash_C31l, nt_hash_D31l, nt_hash_E31l, nt_hash_F31l,
    nt_hash_G31l, nt_hash_H31l, nt_hash_I31l, nt_hash_K31l, nt_hash_L31l,
    nt_hash_M31l, nt_hash_N31l, nt_hash_P31l, nt_hash_Q31l, nt_hash_R31l,
    nt_hash_S31l, nt_hash_T31l, nt_hash_V31l, nt_hash_W31l, nt_hash_Y31l
  };

  static constexpr const std::array<std::array<uint64_t, 33>, 20> msTab33r = {
    nt_hash_A33r, nt_hash_C33r, nt_hash_D33r, nt_hash_E33r, nt_hash_F33r,
    nt_hash_G33r, nt_hash_H33r, nt_hash_I33r, nt_hash_K33r, nt_hash_L33r,
    nt_hash_M33r, nt_hash_N33r, nt_hash_P33r, nt_hash_Q33r, nt_hash_R33r,
    nt_hash_S33r, nt_hash_T33r, nt_hash_V33r, nt_hash_W33r, nt_hash_Y33r
  };

  uint64_t msTab31l_33r_or[20] = {};

  public:
    static constexpr const bool possible_false_positive_matches = true;
    explicit NtHashAminoacidsTransformer(size_t qgram_length)
    {
      for(size_t idx = 0;
          idx < sizeof msTab31l_33r_or/sizeof msTab31l_33r_or[0]; ++idx)
      {
        msTab31l_33r_or[idx] =
          (msTab31l[idx][qgram_length % 31] | msTab33r[idx][qgram_length % 33]);
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
        rev_compl_code ^=
          nt_hash_seed_table[static_cast<size_t>(complement_cc)];
      }
      return std::make_pair(fwd_code,rev_compl_code);
    }

    // forward-strand ntHash for sliding k-mers
    [[nodiscard]] uint64_t next_hash_value_get(uint8_t charOut, uint64_t fhVal,
                                               uint8_t charIn) const noexcept
    {
      assert(charIn < uint8_t{4} && charOut < uint8_t{4});
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
#endif // NTHASH_FWD_AMINOACIDS_HPP
