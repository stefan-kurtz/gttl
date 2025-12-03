#ifndef SKA_LSB_RADIX_SORT_HPP
#define SKA_LSB_RADIX_SORT_HPP
#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include <cassert>
#include <climits>
#include <cstring>
#include <algorithm>
#include <utility>
#include <string>
#include <vector>
#include <tuple>
#include <thread>
#ifndef NDEBUG
#include <iostream>
#endif
#include "utilities/constexpr_for.hpp"
#include "utilities/buckets.hpp"
#include "utilities/lsb_radix_sort.hpp"
#include "utilities/span.hpp"
#include "utilities/runtime_class.hpp"

#ifndef RADIX_SORT8_MAX_SIZEOF_UNIT
#define RADIX_SORT8_MAX_SIZEOF_UNIT 32
#endif

/* code up until function countingsort_skarupke_it was adapted from
   https://probablydance.com/2016/12/27/
             i-wrote-a-faster-sorting-algorithm/#more-4987
*/

template<typename Counttype>
struct PartitionInfo
{
  PartitionInfo()
   : count(0)
  {}

  union
  {
    Counttype count;
    Counttype offset;
  };
};

template<typename Counttype,int sizeof_unit,typename RT>
static inline void swap_two_values(uint8_t *array,Counttype a,Counttype b)
{
  static_assert(sizeof(uint64_t) == 8);
  uint64_t *const aptr = reinterpret_cast<uint64_t *>(array + sizeof_unit * a);
  uint64_t *const bptr = reinterpret_cast<uint64_t *>(array + sizeof_unit * b);
  const uint64_t tmp = *aptr;
  *aptr = *bptr;
  *bptr = tmp;
  RT *const aptr_rest  = reinterpret_cast<RT *>(array + 8 + sizeof_unit * a);
  RT *const bptr_rest  = reinterpret_cast<RT *>(array + 8 + sizeof_unit * b);
  RT tmp_rest = *aptr_rest;
  *aptr_rest = *bptr_rest;
  *bptr_rest = tmp_rest;
}

template <typename Counttype,typename basetype,int sizeof_unit>
static inline void swap_unit(basetype *array,Counttype a,Counttype b)
{
  if constexpr (sizeof_unit == 1)
  {
    std::swap(array[a],array[b]);
  } else
  {
    static_assert(sizeof_unit > 1 && sizeof(basetype) == 1);
    if constexpr (sizeof_unit == 9)
    {
      swap_two_values<Counttype,sizeof_unit,uint8_t>(array,a,b);
    } else
    {
      if constexpr (sizeof_unit == 10)
      {
        swap_two_values<Counttype,sizeof_unit,uint16_t>(array,a,b);
      } else
      {
        if constexpr (sizeof_unit == 12)
        {
          swap_two_values<Counttype,sizeof_unit,uint32_t>(array,a,b);
        } else
        {
          if constexpr (sizeof_unit == 16)
          {
            swap_two_values<Counttype,sizeof_unit,uint64_t>(array,a,b);
          } else
          {
            uint8_t tmp[sizeof_unit];
            memcpy(&tmp[0],&array[sizeof_unit * a],sizeof_unit);
            memcpy(&array[sizeof_unit * a],&array[sizeof_unit * b],sizeof_unit);
            memcpy(&array[sizeof_unit * b],&tmp[0],sizeof_unit);
          }
        }
      }
    }
  }
}

/* This is an iterator over the partitions */
template <typename Numpartitionstype,typename F>
static inline Numpartitionstype *custom_std_partition(Numpartitionstype *begin,
                                                      Numpartitionstype *end,
                                                      F &&func)
{
  for (/* Nothing*/ ; /* Nothing*/; ++begin)
  {
    if (begin == end)
    {
      return end;
    }
    if (!func(*begin))
    {
      break;
    }
  }
  Numpartitionstype *it = begin;
  for (++it; it != end; ++it)
  {
    if (func(*it))
    {
      std::swap(*begin, *it);
      ++begin;
    }
  }
  return begin;
}

template <typename Counttype,typename Func>
static inline void unroll_loop_four_times(Counttype begin,
                                          size_t iteration_count,
                                          Func &&to_call)
{
  size_t loop_count = iteration_count / 4;
  const size_t remainder_count = iteration_count - loop_count * 4;
  for (/* Nothing */; loop_count > 0; --loop_count)
  {
    to_call(begin);
    ++begin;
    to_call(begin);
    ++begin;
    to_call(begin);
    ++begin;
    to_call(begin);
    ++begin;
  }
  if (remainder_count >= 1)
  {
    to_call(begin);
    ++begin;
  }
  if (remainder_count >= 2)
  {
    to_call(begin);
    ++begin;
  }
  if (remainder_count >= 3)
  {
    to_call(begin);
    ++begin;
  }
}

template <typename Counttype,
          typename Numpartitionstype,
          size_t num_buckets,
          typename basetype,
          int sizeof_unit,
          uint64_t (*functor)(int,const basetype *,size_t idx)>
static Buckets<Counttype> *countingsort_skarupke_it(basetype *array,
                                                    size_t array_len,
                                                    int shift)
{
  PartitionInfo<Counttype> partitions[num_buckets];

  for (size_t idx = 0; idx < array_len; ++idx)
  {
    ++partitions[functor(shift,array,idx)].count;
  }
  size_t count_zero = 0;
  for (size_t idx = 0; idx < num_buckets; ++idx)
  {
    count_zero += (partitions[idx].count == 0);
  }
  if (count_zero == num_buckets - 1) /* all but one bucket is empty */
  {
    return nullptr;
  }
  Buckets<Counttype> *buckets = new Buckets<Counttype>(num_buckets);
  Numpartitionstype remaining_partitions[num_buckets];
  Counttype total = 0;
  size_t num_partitions = 0;
  for (size_t idx = 0; idx < num_buckets; ++idx)
  {
    const Counttype count = partitions[idx].count;
    if (count > 0)
    {
      partitions[idx].offset = total;
      total += count;
      remaining_partitions[num_partitions] = idx;
      ++num_partitions;
    }
    buckets->set(idx,total);
  }
  Counttype *bucket_ends = buckets->reference();
  for (Numpartitionstype *last_remaining = remaining_partitions+num_partitions,
                         *end_partition = remaining_partitions + 1;
       last_remaining > end_partition; /* Nothing */)
  {
    last_remaining = custom_std_partition<Numpartitionstype>(
        remaining_partitions,
        last_remaining,
        [&](Numpartitionstype partition)  /* lambda expression */
        {
          Counttype &begin_offset = partitions[partition].offset;
          Counttype &end_offset = bucket_ends[partition];
          if (begin_offset == end_offset)
          {
            return false;
          }

          unroll_loop_four_times<Counttype>(
              begin_offset,
              end_offset - begin_offset,
              [partitions = partitions, array, &shift](Counttype idx)/*lambda*/
              {
                const uint64_t this_partition
                  = functor(shift,array,static_cast<size_t>(idx));
                const Counttype offset = partitions[this_partition].offset++;
                swap_unit<Counttype,basetype,sizeof_unit>(array,idx,offset);
              });
          return begin_offset != end_offset;
        });
  }
#ifndef NDEBUG
  for (auto &&bck : *buckets)
  {
    for (Counttype j = std::get<0>(bck) + 1;  j < std::get<1>(bck); j++)
    {
      const uint64_t a = functor(shift,array,static_cast<size_t>(j-1));
      const uint64_t b = functor(shift,array,static_cast<size_t>(j));
      if (a != b)
      {
        std::cerr << "a=" << a << " != " << b << '\n';
        exit(EXIT_FAILURE);
      }
    }
  }
#endif
  assert(buckets != nullptr);
  return buckets;
}

template<typename Counttype>
// NOLINTNEXTLINE(readability-suspicious-call-argument)
static const Buckets<Counttype> *countingsort_skarupke(
                                                 [[maybe_unused]] int num_bits,
                                                 int shift,
                                                 uint64_t *array,
                                                 size_t array_len)
{
  assert(num_bits == first_pass_msb_bits);
  constexpr const size_t num_buckets = size_t(1) << first_pass_msb_bits;
  const Buckets<Counttype> *const buckets  = countingsort_skarupke_it<
                                Counttype,
                                uint8_t,
                                num_buckets,
                                uint64_t,
                                1,
                                radix_key_uint64<first_pass_msb_bits>>(
                               array, array_len, shift);
  return buckets;
}

template<typename Counttype>
static Buckets<Counttype> *countingsort_skarupke_bytes(int sizeof_unit,
                                                       int byte_index,
                                                       uint8_t *array,
                                                       size_t num_units)
{
  using basetype = uint8_t;
  constexpr const size_t num_buckets = size_t(1) << first_pass_msb_bits;
  assert(byte_index < sizeof_unit);
  assert(sizeof_unit >= 2 && sizeof_unit <= RADIX_SORT8_MAX_SIZEOF_UNIT);
  Buckets<Counttype> *buckets = nullptr;
  constexpr_for<2,RADIX_SORT8_MAX_SIZEOF_UNIT+1,1>([&](auto const_expr_idx)
  {
    if (sizeof_unit == const_expr_idx)
    {
      if constexpr (const_expr_idx < 9)
      {
        buckets = countingsort_skarupke_it<Counttype,uint8_t,num_buckets,
                                           basetype,const_expr_idx,
                                           radix_key_uint8<const_expr_idx>>
                                          (array,num_units,byte_index);
      } else
      {
        buckets = countingsort_skarupke_it<Counttype,uint16_t,num_buckets,
                                           basetype,const_expr_idx,
                                           radix_key_uint8<const_expr_idx>>
                                          (array,num_units,byte_index);
      }
    }
  });
  return buckets;
}

template<typename Counttype,typename basetype,int sizeof_unit>
static std::pair<const Buckets<Counttype> *,int> radixsort_ska_template(
                                                   int num_sort_bits,
                                                   int byte_index,
                                                   basetype *array,
                                                   size_t num_units,
                                                   /* for sizeof_unit>1 */
                                                   [[maybe_unused]]
                                                   bool reversed_byte_order)
{
  int bits_already_sorted = byte_index * first_pass_msb_bits;
  const Buckets<Counttype> *buckets = nullptr;
  while (true)
  {
    bits_already_sorted += first_pass_msb_bits;
    if constexpr (sizeof_unit == 1)
    {
      static_assert(sizeof(basetype) == 8);
      // NOLINTNEXTLINE(readability-suspicious-call-argument)
      buckets = countingsort_skarupke<Counttype>(first_pass_msb_bits,
                                                 64 - bits_already_sorted,
                                                 array,
                                                 num_units);
    } else
    {
      assert(num_sort_bits >= 8);
      static_assert(sizeof(basetype) == 1);
      buckets = countingsort_skarupke_bytes<Counttype>
                                           (sizeof_unit,
                                            real_byte_index(reversed_byte_order,
                                                            byte_index),
                                            array,
                                            num_units);
      byte_index++;
    }
    if (buckets != nullptr || bits_already_sorted >= num_sort_bits)
    {
      break;
    }
  }
  return {buckets,bits_already_sorted};
}

template<typename Counttype,typename basetype,int sizeof_unit,class SorterClass>
static const Buckets<Counttype> *radixsort_ska_then_other_generic(
                                                   SorterClass *sorter_instance,
                                                   int num_sort_bits,
                                                   basetype *array,
                                                   size_t num_units)
{
  if (num_units < 2)
  {
    return nullptr;
  }
  auto result = radixsort_ska_template<Counttype,basetype,sizeof_unit>
                                      (num_sort_bits,
                                       0,
                                       array,
                                       num_units,
                                       /* for sizeof_unit>1 */
                                       sorter_instance
                                         ->reversed_byte_order_get());
  const int bits_already_sorted = std::get<1>(result);
  const Buckets<Counttype> *const buckets = std::get<0>(result);
  if (buckets == nullptr || bits_already_sorted >= num_sort_bits)
  {
    return buckets;
  }
  const Counttype maximum_bucket_width = buckets->maximum_width_get();
  sorter_instance->setup(maximum_bucket_width);
  for (auto &&bck : *buckets)
  {
    const Counttype bucket_start = std::get<0>(bck);
    const Counttype bucket_width = std::get<1>(bck) - bucket_start;
    if (bucket_width > 1)
    {
      sorter_instance->sort(array,bucket_start,bucket_width,
                            bits_already_sorted);
    }
  }
  return buckets;
}

class LSBuint64Sorter
{
  private:
    uint64_t *uint64_buffer;
    size_t buffer_width;
    const int num_sort_bits;
  public:
  LSBuint64Sorter(int _num_sort_bits)
    : uint64_buffer(nullptr)
    , buffer_width(0)
    , num_sort_bits(_num_sort_bits)
  { }
  ~LSBuint64Sorter(void)
  {
    delete[] uint64_buffer;
  }
  void setup(size_t _bucket_width)
  {
    buffer_width = _bucket_width;
    uint64_buffer = new uint64_t [_bucket_width];
  }
  [[nodiscard]] bool reversed_byte_order_get(void) const noexcept
  {
    return false; /* Not used */
  }
  void sort(uint64_t *array,size_t bucket_start,size_t bucket_width,
            int bits_already_sorted)
  {
    assert(bucket_width <= buffer_width);
    lsb_radix_sort(array + bucket_start,
                   uint64_buffer,
                   bucket_width,
                   bits_already_sorted,
                   num_sort_bits - bits_already_sorted);
  }
};

template<int sizeof_unit>
class LSBbytesSorter
{
  private:
    uint8_t *uint8_buffer;
    const int num_sort_bits;
    const bool reversed_byte_order;
  public:
  LSBbytesSorter(int _num_sort_bits,
                 bool _reversed_byte_order)
    : uint8_buffer(nullptr)
    , num_sort_bits(_num_sort_bits)
    , reversed_byte_order(_reversed_byte_order)
   { }
  ~LSBbytesSorter(void)
  {
    delete[] uint8_buffer;
  }
  void setup(size_t maximum_bucket_width)
  {
    uint8_buffer = new uint8_t [sizeof_unit * maximum_bucket_width];
  }
  [[nodiscard]] bool reversed_byte_order_get(void) const noexcept
  {
    return reversed_byte_order;
  }
  void sort(uint8_t *array,size_t bucket_start,size_t bucket_width,
            int bits_already_sorted)
  {
    lsb_radix_sort<sizeof_unit>(array + sizeof_unit * bucket_start,
                                uint8_buffer,
                                bucket_width,
                                bits_already_sorted,
                                num_sort_bits - bits_already_sorted,
                                reversed_byte_order);
  }
};

template<typename Counttype>
static const Buckets<Counttype> *ska_lsb_radix_sort(int num_sort_bits,
                                              uint64_t *array,
                                              size_t array_len)
{
  LSBuint64Sorter lsb_uint64_sorter(num_sort_bits);
  return radixsort_ska_then_other_generic<Counttype,uint64_t,1,LSBuint64Sorter>
                                          (&lsb_uint64_sorter,
                                           num_sort_bits,
                                           array,
                                           array_len);
}

template<typename Counttype>
static Buckets<Counttype> *ska_lsb_radix_sort(int sizeof_unit,
                                              int num_sort_bits,
                                              uint8_t *array,
                                              size_t num_units,
                                              bool reversed_byte_order)
{
  assert(sizeof_unit >= 2 && sizeof_unit <= RADIX_SORT8_MAX_SIZEOF_UNIT);
  assert(sizeof_unit * CHAR_BIT >= num_sort_bits);
  Buckets<Counttype> *buckets = nullptr;
  constexpr_for<2,RADIX_SORT8_MAX_SIZEOF_UNIT+1,1>([&](auto const_expr_idx)
  {
    if (sizeof_unit == const_expr_idx)
    {
      LSBbytesSorter<const_expr_idx>
         lsb_bytes_sorter(num_sort_bits,reversed_byte_order);
      buckets = radixsort_ska_then_other_generic<Counttype,uint8_t,
                                                 const_expr_idx,
                                                 LSBbytesSorter<const_expr_idx>>
                                                 (lsb_bytes_sorter,
                                                  num_sort_bits,
                                                  array,
                                                  num_units);
    }
  });
  return buckets;
}

#undef SHOW_NON_EMPTY_BUCKETS
#ifdef SHOW_NON_EMPTY_BUCKETS
static void show_non_empty_buckets(const Buckets<size_t> *buckets,
                                   int byte_index,size_t num_units)
{
  size_t idx = 0;
  for (auto it = buckets->begin(); it != buckets->end(); ++it)
  {
    size_t bucket_width = std::get<1>(*it) - std::get<0>(*it);
    if (bucket_width > 0)
    {
      for (int i = 0; i < byte_index; i++)
      {
        printf("%c",'\t');
      }
      printf("bucket\t%03lu\t%zu\n",idx,bucket_width);
    }
    idx++;
  }
  for (int i = 0; i < byte_index; i++)
  {
    printf("%c",'\t');
  }
  printf("maximum width of buckets\t%zu\tnum_unit\t%zu\n",
         buckets->maximum_width_get(),num_units);
}
#endif

template<typename Counttype,typename basetype,int sizeof_unit,class SorterClass>
static void ska_large_lsb_small_radix_sort_generic(SorterClass *sorter_instance,
                                                   int num_sort_bits,
                                                   basetype *array,
                                                   size_t num_units)
{
  if (num_units < 2)
  {
    return;
  }
  const size_t skarupke_threshold = num_units/10;
  const int num_sort_bytes = (num_sort_bits+CHAR_BIT-1)/CHAR_BIT;
  using StackStruct = struct { size_t offset;
                               size_t num_units;
                               int byte_index;
                             };

  std::vector<StackStruct> stack{{size_t(0), num_units, 0}};
  sorter_instance->setup(skarupke_threshold);
  while (not stack.empty())
  {
    const StackStruct current = stack.back();
    stack.pop_back();

    const Buckets<Counttype> *buckets;
    if constexpr (sizeof_unit == 1)
    {
      static_assert(sizeof(basetype) == 8);
      buckets = countingsort_skarupke<Counttype>(first_pass_msb_bits,
                                                 56 -
                                                 CHAR_BIT * current.byte_index,
                                                 array + current.offset,
                                                 current.num_units);
    } else
    {
      assert(num_sort_bits >= 8);
      const int this_byte_index = real_byte_index(sorter_instance
                                                    ->reversed_byte_order_get(),
                                                  current.byte_index);
      static_assert(sizeof(basetype) == 1);
      buckets = countingsort_skarupke_bytes<Counttype>
                                           (sizeof_unit,
                                            this_byte_index,
                                            array +
                                              current.offset * sizeof_unit,
                                            current.num_units);
    }
    if (current.byte_index + 1 < num_sort_bytes)
    {
      if (buckets == nullptr) /* all elements in one bucket */
      {
        stack.push_back({current.offset,
                         current.num_units,
                         current.byte_index + 1});
      } else
      {
#ifdef SHOW_NON_EMPTY_BUCKETS
        show_non_empty_buckets(buckets,current.byte_index,current.num_units);
#endif
        for (auto &bck : *buckets)
        {
          const Counttype bucket_start = std::get<0>(bck);
          const Counttype bucket_width = std::get<1>(bck) - bucket_start;
          if (bucket_width > skarupke_threshold)
          {
            stack.push_back({current.offset + bucket_start,
                             static_cast<size_t>(bucket_width),
                             current.byte_index + 1});
          } else
          {
            if (bucket_width > 1)
            {
              const int bits_already_sorted = CHAR_BIT * current.byte_index;
              sorter_instance->sort(array,
                                    current.offset + bucket_start,
                                    bucket_width,
                                    bits_already_sorted);
            }
          }
        }
      }
    }
    delete buckets;
  }
}

class PartInfoTab
{
  private:
  std::vector<size_t> end_indexes;
  public:
  PartInfoTab(size_t num_elements,size_t num_parts)
  {
    if (num_parts > num_elements)
    {
      end_indexes.push_back(num_elements);
    } else
    {
      const size_t avg_width = (num_elements + num_parts - 1)/num_parts;
      for (size_t p = 0; p < num_parts; p++)
      {
        end_indexes.push_back(std::min((p+1) * avg_width,num_elements));
      }
    }
  }
  [[nodiscard]] size_t size(void) const noexcept { return end_indexes.size(); }
  std::pair<size_t,size_t> operator [](size_t idx) const noexcept
  {
    return std::make_pair(idx == 0 ? 0 : end_indexes[idx-1],
                          idx == 0 ? end_indexes[0]
                                   : (end_indexes[idx] - end_indexes[idx-1]));
  }
};

#ifndef NDEBUG
template<class It>
static inline void check_order(It begin, It end)
{
  auto previous = *begin;
  for(auto it = begin + 1; it != end; ++it)
  {
    if (previous > *it)
    {
      assert(false);
    }
    previous = *it;
  }
}
#endif

static inline void ska_large_lsb_small_radix_sort(int num_sort_bits,
                                                  uint64_t *array,
                                                  size_t num_units,
                                                  size_t num_threads = 1)
{
  using Counttype = size_t;
  if (num_threads == 1)
  {
    LSBuint64Sorter lsb_uint64_sorter(num_sort_bits);
    ska_large_lsb_small_radix_sort_generic<Counttype,uint64_t,1,LSBuint64Sorter>
                                          (&lsb_uint64_sorter,
                                           num_sort_bits,
                                           array,
                                           num_units);
  } else
  {
    PartInfoTab part_info_tab(num_units,num_threads);
    std::vector<std::thread> threads;
    threads.reserve(part_info_tab.size());
    for (size_t thd_num = 0; thd_num < part_info_tab.size(); thd_num++)
    {
      threads.emplace_back([&part_info_tab, array,
                            num_sort_bits, thd_num]
      {
        size_t left;
        size_t width;
        std::tie(left,width) = part_info_tab[thd_num];
        LSBuint64Sorter lsb_uint64_sorter(num_sort_bits);
        ska_large_lsb_small_radix_sort_generic<Counttype,uint64_t,1,
                                               LSBuint64Sorter>
                                              (&lsb_uint64_sorter,
                                               num_sort_bits,
                                               array + left,
                                               width);
      });
    }
    for (auto &th : threads)
    {
      th.join();
    }
    RunTimeClass rt_multi_way_merge{};
    Span<uint64_t> span(array,num_units);
    for (size_t idx = 1; idx < part_info_tab.size(); idx++)
    {
      size_t left;
      size_t width;
      std::tie(left,width) = part_info_tab[idx];
#ifndef NDEBUG
      check_order(span.begin(),span.begin() + left);
      check_order(span.begin() + left, span.begin() + left + width);
#endif
      std::inplace_merge(span.begin(),
                         span.begin() + left,
                         span.begin() + left + width);
    }
    rt_multi_way_merge.show(std::string("merging ") +
                            std::to_string(part_info_tab.size()) +
                            std::string(" sorted vectors"));
  }
}

static inline void ska_large_lsb_small_radix_sort(int sizeof_unit,
                                                  int num_sort_bits,
                                                  uint8_t *array,
                                                  size_t num_units,
                                                  bool reversed_byte_order)
{
  using Counttype = size_t;
  assert(sizeof_unit >= 2 and
         sizeof_unit <= RADIX_SORT8_MAX_SIZEOF_UNIT and
         sizeof_unit * CHAR_BIT >= num_sort_bits);
  constexpr_for<2,RADIX_SORT8_MAX_SIZEOF_UNIT+1,1>([&](auto const_expr_idx)
  {
    if (sizeof_unit == const_expr_idx)
    {
      LSBbytesSorter<const_expr_idx> lsb_bytes_sorter(num_sort_bits,
                                                      reversed_byte_order);
      ska_large_lsb_small_radix_sort_generic<Counttype,uint8_t,const_expr_idx,
                                             LSBbytesSorter<const_expr_idx>>
                                            (&lsb_bytes_sorter,
                                             num_sort_bits,
                                             array,
                                             num_units);
    }
  });
}
#endif
