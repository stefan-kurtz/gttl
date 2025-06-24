#ifndef NTHASH_ROTATION_TABLES_HPP
#define NTHASH_ROTATION_TABLES_HPP

#include <array>
#include <cstddef>
#include <cstdint>

constexpr uint64_t MASK_33 = (uint64_t{1} << uint64_t{33}) - 1;
constexpr uint64_t MASK_TOP_31 = ~((uint64_t{1} << uint64_t{33}) - 1);

consteval uint64_t rolR33(uint64_t value)
{
  // Using a shorter type here would lead to implicit conversions
  const uint64_t first_bit = value >> uint64_t{32};
  return ((value << uint64_t{1}) | first_bit) & MASK_33;
}

consteval uint64_t rolL31(uint64_t value)
{
  const uint64_t first_bit = value >> uint64_t{63};
  const uint64_t moved_bit = first_bit << uint64_t{33};
  return ((value << uint64_t{1}) | moved_bit) & MASK_TOP_31;
}

consteval std::array<uint64_t, 33> generate_33r_table(uint64_t seed)
{
  std::array<uint64_t, 33> table{};

  uint64_t last_33_bits = seed & MASK_33;
  for(size_t i = 0; i < 33; ++i)
  {
    table[i] = last_33_bits;
    last_33_bits = rolR33(last_33_bits);
  }
  return table;
}

consteval std::array<uint64_t, 31> generate_31l_table(uint64_t seed)
{
  std::array<uint64_t, 31> table{};

  uint64_t first_31_bits = seed & MASK_TOP_31;
  for(size_t i = 0; i < 31; ++i)
  {
    table[i] = first_31_bits;
    first_31_bits = rolL31(first_31_bits);
  }
  return table;
}


#endif // NTHASH_ROTATION_TABLES_HPP
