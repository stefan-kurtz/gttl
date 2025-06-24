#ifndef NTHASH_FWD_AMINOACIDS_HPP
#define NTHASH_FWD_AMINOACIDS_HPP

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <utility>
#include "sequences/nthash_rotation_tables.hpp"

template <size_t N, size_t... I>
consteval static auto make31l_table_inner(const std::array<uint64_t, N>& input,
                                          std::index_sequence<I...> /**/)
{
  return std::array{ generate_31l_table(input[I])... };
}

template <size_t N>
consteval static auto make31l_table(const std::array<uint64_t, N>& input)
{
  return make31l_table_inner(input, std::make_index_sequence<N>{});
}

template <size_t N, size_t... I>
consteval static auto make33r_table_inner(const std::array<uint64_t, N>& input,
                                          std::index_sequence<I...> /**/)
{
  return std::array { generate_33r_table(input[I])... };
}

template <size_t N>
consteval static auto make33r_table(const std::array<uint64_t, N>& input)
{
  return make33r_table_inner(input, std::make_index_sequence<N>{});
}

class NtHashAminoacidsTransformer
{
  // rotate "v" to the left 1 position
  [[nodiscard]] static uint64_t rotate_left_1(const uint64_t v) noexcept
  {
    return (v << 1) | (v >> 63);
  }

  // swap bit 0 with bit 33 in "v"
  [[nodiscard]] static uint64_t swapbits033(const uint64_t v) noexcept
  {
    const uint64_t x = (v ^ (v >> 33)) & uint64_t{1};
    return v ^ (x | (x << 33));
  }

  static constexpr const std::array<uint64_t, 20> nt_hash_seed_table = {
    0xE6CC5634828DE226, 0xC3916818AB69EAEF, 0x708B967AB82D90DA,
    0x97EE132403084238, 0x54B5F4F96AA2B4C6, 0xA16E4815AEB82F9C,
    0xD45994107632BA9E, 0x52CA8A666024FEB4, 0x4B7BE4E0445E0FB4,
    0x10936D371508CFBB, 0x9845685EF8E90AC2, 0xBC57ED4598240F30,
    0x4ADD75731C811286, 0xF4EADE76635CE527, 0x725F7D219B2873BC,
    0x2A90195B4166C058, 0xEABF48D47DA98AED, 0x7BE88554E92B4BF9,
    0x211DFFCF403D167C, 0x8454A6A3B239AE9E
  };

  static constexpr const auto msTab31l = make31l_table(nt_hash_seed_table);
  static constexpr const auto msTab33r = make33r_table(nt_hash_seed_table);

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
        assert(t_sequence[idx] < uint8_t(20));
        hVal ^= nt_hash_seed_table[t_sequence[idx]];
      }
      return hVal;
    }

    // forward-strand ntHash for sliding k-mers
    [[nodiscard]] uint64_t next_hash_value_get(uint8_t charOut, uint64_t fhVal,
                                               uint8_t charIn) const noexcept
    {
      assert(charIn < uint8_t(20) && charOut < uint8_t(20));
      uint64_t hVal = rotate_left_1(fhVal);
      hVal = swapbits033(hVal);
      hVal ^= nt_hash_seed_table[charIn];
      hVal ^= msTab31l_33r_or[charOut];
      return hVal;
    }
};
#endif // NTHASH_FWD_AMINOACIDS_HPP
