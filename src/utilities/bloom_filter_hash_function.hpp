/*
** Developed by Henning Lindemann
 */

#ifndef BLOOM_FILTER_HASH_FUNCTION_HPP
#define BLOOM_FILTER_HASH_FUNCTION_HPP

#include <cstddef>
#include <cstdint>

static inline uint64_t hash_function(uint64_t value, size_t seed)
{
  static constexpr const uint64_t FNV_prime = 0x00000100000001b3;
  uint64_t hash_value = 0xcbf29ce484222325;

  hash_value ^= static_cast<uint64_t>(seed);
  hash_value *= FNV_prime;

  for (uint_fast8_t idx = 0; idx < 8; idx++)
  {
    hash_value ^= ((value >> (idx * 8U)) & 0xFFU);
    hash_value *= FNV_prime;
  }
  return hash_value;
}


#endif // BLOOM_FILTER_HASH_FUNCTION_HPP
