#ifndef LSB_RADIX_SORT_HPP
#define LSB_RADIX_SORT_HPP
#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include <cassert>
#include <climits>
#include <cstring>
#include <algorithm>
#include "utilities/mathsupport.hpp"
#include "utilities/bits2split.hpp"

static inline int real_byte_index(bool reversed_byte_order,int byte_index)
{
  if (reversed_byte_order && byte_index < 8)
  {
    return 8 - 1 - byte_index;
  }
  return byte_index;
}

template<typename basetype,int num_bits,int sizeof_unit,
         uint64_t (*functor)(int,const basetype *,size_t)>
static bool lsb_radix_sort_single_pass(
                                int functor_param1, /* shift or byte_index */
                                basetype *dest,
                                const basetype *src,
                                size_t array_len)
{
  constexpr const size_t num_buckets = size_t(1) << num_bits;
  size_t count_zero = 0;
  size_t count[num_buckets] = {0};

  assert(array_len > 1);
  for (size_t idx = 0; idx < array_len; ++idx)
  {
    ++count[functor(functor_param1,src,idx)];
  }
  for (size_t idx = 0; idx < num_buckets; ++idx)
  {
    count_zero += (count[idx] == 0);
  }
  if (count_zero == num_buckets - 1)
  {
    return false;
  }
  for (size_t idx = 1; idx < num_buckets; ++idx)
  {
    count[idx] += count[idx-1];
  }
  size_t idx = array_len;
  do
  {
    idx--;
    const uint64_t key = functor(functor_param1,src,idx);
    const size_t d = --count[key];
    if constexpr (sizeof_unit == 1)
    {
      dest[d] = src[idx];
    } else
    {
      memcpy(dest + sizeof_unit * d,src + sizeof_unit * idx,sizeof_unit);
    }
  } while (idx > 0);
  return true;
}

template<int num_bits>
static inline uint64_t radix_key_uint64(int shift,const uint64_t *array,
                                        size_t idx)
{
  constexpr const uint64_t mask = gttl_bits2maxvalue<uint64_t,num_bits>();
  return (array[idx] >> shift) & mask;
}

template<int sizeof_unit>
static inline uint64_t radix_key_uint8(int byte_index,const uint8_t *array,
                                       size_t idx)
{
  return static_cast<uint64_t>(array[sizeof_unit * idx + byte_index]);
}

static constexpr const int first_pass_msb_bits = 8;

template<typename basetype>
using SinglePassSorter = bool (*) (int, basetype *,const basetype *, size_t);

#define LSB_RADIX_SORT_SINGLE_PASS(NUM_BITS)\
        lsb_radix_sort_single_pass<uint64_t,NUM_BITS,1,\
                                   radix_key_uint64<NUM_BITS>>

static inline void lsb_radix_sort(uint64_t *array,
                                  uint64_t *buffer,
                                  size_t array_len,
                                  int bits_already_sorted,
                                  int remaining_bits)
{
  static constexpr const SinglePassSorter<uint64_t> sorter_table_uint64[] =
  {
    LSB_RADIX_SORT_SINGLE_PASS(1),
    LSB_RADIX_SORT_SINGLE_PASS(2),
    LSB_RADIX_SORT_SINGLE_PASS(3),
    LSB_RADIX_SORT_SINGLE_PASS(4),
    LSB_RADIX_SORT_SINGLE_PASS(5),
    LSB_RADIX_SORT_SINGLE_PASS(6),
    LSB_RADIX_SORT_SINGLE_PASS(7),
    LSB_RADIX_SORT_SINGLE_PASS(8),
    LSB_RADIX_SORT_SINGLE_PASS(9)
  };
  size_t num_ranges;
  const int *const bit_groups  = bit_counts_get(&num_ranges, remaining_bits);
  constexpr const int max_bits = 64;
  assert(bits_already_sorted + remaining_bits <= max_bits);
  int shift = max_bits - (bits_already_sorted + remaining_bits);
  // NOLINTBEGIN(misc-const-correctness) since linters do not identify std::swap
  uint64_t *src_ptr = array;
  uint64_t *dest_ptr = buffer;
  // NOLINTEND(misc-const-correctness)
  for (size_t idx = 0; idx < num_ranges; idx++)
  {
    const int bits = bit_groups[idx];
    const SinglePassSorter<uint64_t> func = sorter_table_uint64[bits - 1];
    const bool was_permuted = func(shift,dest_ptr,src_ptr,array_len);
    if (was_permuted)
    {
      std::swap(src_ptr,dest_ptr);
    }
    shift += bits;
  }
  assert(shift == max_bits - bits_already_sorted);
  if (src_ptr != array)
  {
    memcpy(array,buffer,array_len * sizeof *array);
  }
}

template<int sizeof_unit>
static void lsb_radix_sort(uint8_t *array,
                           uint8_t *buffer,
                           size_t array_len,
                           int bits_already_sorted,
                           int remaining_bits,
                           bool reversed_byte_order)
{
  const int last_byte_index = remaining_bits / CHAR_BIT +
                              ((remaining_bits % CHAR_BIT) != 0 ? 1 : 0);
  const SinglePassSorter<uint8_t> func = lsb_radix_sort_single_pass<
                               uint8_t,
                               first_pass_msb_bits,
                               sizeof_unit,
                               radix_key_uint8<sizeof_unit>>;
  // NOLINTBEGIN(misc-const-correctness)
  uint8_t *src_ptr = array;
  uint8_t *dest_ptr = buffer;
  // NOLINTEND(misc-const-correctness)
  assert(bits_already_sorted % 8 == 0);
  for (int byte_index = last_byte_index;
       byte_index >= bits_already_sorted/8; byte_index--)
  {
    const bool was_permuted = func(real_byte_index(reversed_byte_order,
                                                   byte_index),
                                   dest_ptr,src_ptr,array_len);
    if (was_permuted)
    {
      std::swap(src_ptr,dest_ptr);
    }
  }
  if (src_ptr != array)
  {
    memcpy(array,buffer,sizeof_unit * array_len * sizeof *array);
  }
}
#endif
