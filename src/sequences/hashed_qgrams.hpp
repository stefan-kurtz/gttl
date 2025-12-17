#ifndef HASHED_QGRAMS_HPP
#define HASHED_QGRAMS_HPP
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <stdexcept>
#include <iostream>
#include <string>
#include <cstring> /* for cache input/output */
#include <fstream> /* for cache input/output */
#include <utility>
#include <tuple>
#include <vector>
#include <algorithm>
#include <cmath>
#include <format>
#include <deque>

#include "utilities/bitpacker.hpp"
#include "utilities/buckets.hpp"
#include "utilities/ska_lsb_radix_sort.hpp"
#include "utilities/mathsupport.hpp"
#include "utilities/runtime_class.hpp"
#include "utilities/bytes_unit.hpp"
#include "utilities/is_big_endian.hpp"
#include "threading/thread_pool_var.hpp"
#include "sequences/char_range.hpp"
#include "sequences/char_finder.hpp"
#include "sequences/gttl_multiseq.hpp"

template<int sizeof_unit>
using HashedQgramVector = std::vector<BytesUnit<sizeof_unit,3>>;

/* Header for cache of reference index,
   so the load-method can identify it */
struct HashedQgramsCacheHeader
{
  char magic[8];
  uint32_t version;
  uint8_t sizeof_unit;
  uint8_t handle_both_strands;
  uint8_t possible_false_positive;
  uint8_t at_constant_distance;
  uint8_t has_wildcards;
  uint8_t reserved[3];
  uint64_t qgram_length;
  int64_t hashbits;
  uint64_t count_all_qgrams;
  uint64_t ref_window_size;
  uint64_t sequences_number_bits;
  uint64_t sequences_length_bits;
  uint64_t hashed_qram_vector_size;
};

static constexpr const char hashed_qgrams_cache_magic[8]
  = {'H','Q','G','I','D','X','0','1'};
static constexpr const uint32_t hashed_qgrams_cache_version = 1;

/* validation is slow, so make sure that #undef is used */
#undef VALIDATE_MINIMIZER
#ifdef VALIDATE_MINIMIZER
#include <set>

template<int sizeof_unit>
static void validate_minimizers(
                size_t window_size,
                const HashedQgramVector<sizeof_unit> &all_hashed_qgrams,
                const HashedQgramVector<sizeof_unit> *minimizer_vector,
                const GttlBitPacker<sizeof_unit,3> &hashed_qgram_packer,
                size_t search_start)
{
  std::set<uint64_t> minimizer_set{};
  for (auto it = minimizer_vector->begin() + search_start;
       it != minimizer_vector->end(); it++)
  {
    const uint64_t this_hash = (*it).template decode_at<0>(hashed_qgram_packer);
    (void) minimizer_set.insert(this_hash);
  }
  if (window_size > all_hashed_qgrams.size())
  {
    return;
  }
  for (size_t widx = 0; widx < all_hashed_qgrams.size() - window_size + 1;
       widx++)
  {
    uint64_t min_hash
      = all_hashed_qgrams[widx].template decode_at<0>(hashed_qgram_packer);
    for (size_t j = widx + 1; j < widx + window_size; j++)
    {
      const uint64_t this_hash
        = all_hashed_qgrams[j].template decode_at<0>(hashed_qgram_packer);
      if (this_hash < min_hash)
      {
        min_hash = this_hash;
      }
    }
    if (minimizer_set.count(min_hash) == 0)
    {
      throw std::runtime_error(std::format("cannot find minimizer {}",
                                           min_hash));
    }
  }
}
#endif

static constexpr const char_finder::NucleotideFinder unw_nucleotide_finder{};
using NucleotideRanger = GttlCharRange<char_finder::NucleotideFinder,
                                       unw_nucleotide_finder,
                                       true, false>;

template<int sizeof_unit,class HashIterator>
static std::pair<size_t,bool> append_minimizers(
                                       size_t qgram_length,
                                       size_t window_size,
                                       uint64_t hash_mask,
                                       const GttlBitPacker<sizeof_unit,3>
                                         &hashed_qgram_packer,
                                       HashedQgramVector<sizeof_unit>
                                         *minimizer_vector,
                                       const char *sequence,
                                       size_t seqlen,
                                       size_t seqnum)
{
  if constexpr (HashIterator::handle_both_strands)
  {
    if (seqnum % 2 == 1) /* for processing the reverse complement
                            we skip every second sequence, as the
                            minimizers come from the original sequence */
    {
      return std::make_pair(0, false);
    }
  }
  const size_t minseqlen_to_process = window_size + qgram_length - 1;
  std::deque<BytesUnit<sizeof_unit,3>> window_deque{};
  size_t count_all_qgrams = 0;
  bool has_wildcards = false;
  const NucleotideRanger ranger(sequence, seqlen);
  HashedQgramVector<sizeof_unit>  //NOLINT(misc-const-correctness)
    palindromic_vector{};
  for (auto const &&range : ranger)
  {
    const size_t this_length = std::get<1>(range);
    has_wildcards = has_wildcards || this_length < seqlen;
    if (this_length < minseqlen_to_process)
    {
      continue;
    }
    size_t seqpos = std::get<0>(range);
    const char *const seqptr = sequence + seqpos;
    HashIterator qgiter(qgram_length, seqptr, this_length);
    count_all_qgrams += (this_length - qgram_length - 1);

#ifdef VALIDATE_MINIMIZER
    HashedQgramVector<sizeof_unit> all_hashed_qgrams;
    size_t search_start = 0;
#endif
    bool front_was_moved = false;
    for (auto const &&code_pair : qgiter)
    {
      uint64_t this_hash = // NOLINT(misc-const-correctness)
        std::get<0>(code_pair) & hash_mask;
      size_t stored_seqnum;
      size_t stored_seqpos;
      if constexpr (HashIterator::handle_both_strands)
      {
        const uint64_t rc_hash = std::get<1>(code_pair) & hash_mask;
        if (rc_hash < this_hash)
        {
          this_hash = rc_hash;
          stored_seqnum = seqnum + 1;
          assert(seqlen >= seqpos + qgram_length);
          stored_seqpos = seqlen - (seqpos + qgram_length);
        } else
        {
          if (rc_hash == this_hash)
          {
            /* also add the coordinates for the reverse complement
               as it would otherwise be neglected */
            const BytesUnit<sizeof_unit, 3> palindromic_hashed_qgram_rc(
                                         hashed_qgram_packer,
                                         {this_hash,
                                          static_cast<uint64_t>(seqnum + 1),
                                          static_cast<uint64_t>(
                                                             seqlen
                                                             - (seqpos
                                                             + qgram_length))});
            palindromic_vector.emplace_back(palindromic_hashed_qgram_rc);
            const BytesUnit<sizeof_unit, 3> palindromic_hashed_qgram_fwd(
                                         hashed_qgram_packer,
                                         {this_hash,
                                          static_cast<uint64_t>(seqnum),
                                          static_cast<uint64_t>(seqpos)});
            palindromic_vector.emplace_back(palindromic_hashed_qgram_fwd);
          }
          stored_seqnum = seqnum;
          stored_seqpos = seqpos;
        }
      } else
      {
        stored_seqnum = seqnum;
        stored_seqpos = seqpos;
      }
      while (true)
      {
        if (window_deque.empty())
        {
          front_was_moved = false; /* as next element becomes new front
                                      which has not been moved before */
          break;
        }
        if (window_deque.back().template decode_at<0>(hashed_qgram_packer)
            <= this_hash)
        {
          break;
        }
        window_deque.pop_back();
      }
      const BytesUnit<sizeof_unit, 3> current_hashed_qgram(
                                   hashed_qgram_packer,
                                   {this_hash,
                                    static_cast<uint64_t>(stored_seqnum),
                                    static_cast<uint64_t>(stored_seqpos)});
#ifdef VALIDATE_MINIMIZER
      all_hashed_qgrams.emplace_back(current_hashed_qgram);
#endif
      window_deque.emplace_back(current_hashed_qgram);

      // After the minimizer of the first window was found
      // At this point we are only looking if
      // the current minimizer drops out of the window
      if (seqpos >= window_size)
      {
        uint64_t min_pos_in_window = //NOLINT(misc-const-correctness)
          window_deque.front().template decode_at<2>
          (hashed_qgram_packer);
        if constexpr (HashIterator::handle_both_strands)
        {
          if (window_deque.front().template decode_at<1>(hashed_qgram_packer) ==
              seqnum + 1)
          {
            min_pos_in_window = seqlen - (min_pos_in_window + qgram_length);
          }
        }
        // check if current minimizer drops out of the window
        if (min_pos_in_window <= seqpos - window_size)
        {
          window_deque.pop_front();
          front_was_moved = false; /* queue is still not empty and
                                      now front element was not moved yet */
        }
        assert(!window_deque.empty()); /* because we added
                                          current_hashed_qgram */
        /* check if new minimizer was found (i.e. front element which was
           not moved yet */
        if (!front_was_moved)
        {
          minimizer_vector->emplace_back(window_deque.front());
          front_was_moved = true; /* we moved front element and do not
                                     want to do it again */
        }
      } else
      {
        // add minimizer of first window
        if (seqpos == window_size - 1)
        {
          minimizer_vector->emplace_back(window_deque.front());
          front_was_moved = true;
        }
      }
      seqpos++;
    }
#ifdef VALIDATE_MINIMIZER
    validate_minimizers<sizeof_unit>(window_size,
                                     all_hashed_qgrams,
                                     minimizer_vector,
                                     hashed_qgram_packer,
                                     search_start);
#endif
    window_deque.clear();
  }
  if constexpr (HashIterator::handle_both_strands)
  {
    for (auto && pq : palindromic_vector)
    {
      minimizer_vector->emplace_back(pq);
    }
  }
  return std::make_pair(count_all_qgrams, has_wildcards);
}

template<int sizeof_unit,class HashIterator>
static std::pair<size_t,bool> append_constant_distance_hashed_qgrams(
                                       size_t qgram_length,
                                       size_t window_size,
                                       uint64_t hash_mask,
                                       const GttlBitPacker<sizeof_unit,3>
                                         &hashed_qgram_packer,
                                       HashedQgramVector<sizeof_unit>
                                         *hashed_qgrams_vector,
                                       const char *sequence,
                                       size_t seqlen,
                                       size_t seqnum)
{
  if constexpr (HashIterator::handle_both_strands)
  {
    if (seqnum % 2 == 1) /* for processing the reverse complement
                            we skip every second sequence, as the
                            minimizers come from the original sequence */
    {
      return std::make_pair(0, false);
    }
  }
  const size_t minseqlen_to_process = window_size + qgram_length - 1;
  size_t count_all_qgrams = 0;
  bool has_wildcards = false;
  const NucleotideRanger ranger(sequence, seqlen);
  for (auto const &&range : ranger)
  {
    const size_t this_length = std::get<1>(range);
    has_wildcards = has_wildcards || this_length < seqlen;
    if (this_length < minseqlen_to_process)
    {
      continue;
    }
    size_t seqpos = std::get<0>(range);
    const char *const seqptr = sequence + seqpos;
    HashIterator qgiter(qgram_length, seqptr, this_length);
    count_all_qgrams += (this_length - qgram_length - 1);
    size_t steps = 0;
    for (auto const &&code_pair : qgiter)
    {
      if (steps == 0)
      {
        uint64_t this_hash = //NOLINT(misc-const-correctness)
          std::get<0>(code_pair) & hash_mask;
        size_t stored_seqnum;
        size_t stored_seqpos;
        if constexpr (HashIterator::handle_both_strands)
        {
          const uint64_t rc_hash = std::get<1>(code_pair) & hash_mask;
          if (rc_hash < this_hash)
          {
            this_hash = rc_hash;
            stored_seqnum = seqnum + 1;
            assert(seqlen >= seqpos + qgram_length);
            stored_seqpos = seqlen - (seqpos + qgram_length);
          } else
          {
            stored_seqnum = seqnum;
            stored_seqpos = seqpos;
          }
        } else
        {
          stored_seqnum = seqnum;
          stored_seqpos = seqpos;
        }
        const BytesUnit<sizeof_unit, 3> current_hashed_qgram(
                                     hashed_qgram_packer,
                                     {this_hash,
                                      static_cast<uint64_t>(stored_seqnum),
                                      static_cast<uint64_t>(stored_seqpos)});
        hashed_qgrams_vector->emplace_back(current_hashed_qgram);
        steps = window_size;
      }
      assert(steps > 0);
      steps--;
      seqpos++;
    }
  }
  return std::make_pair(count_all_qgrams, has_wildcards);
}

template<int sizeof_unit>
struct HashedQgramVectorTable
{
  /* common data */
  std::vector<size_t> count_all_qgrams;
  std::vector<std::atomic<bool>> has_wildcards;
  std::vector<HashedQgramVector<sizeof_unit>> table;
  HashedQgramVectorTable(size_t number_of_threads)
    : count_all_qgrams(number_of_threads,0)
    , table(number_of_threads,HashedQgramVector<sizeof_unit>{})
  {
    has_wildcards = std::vector<std::atomic<bool>>(number_of_threads);
    for(auto &b : has_wildcards) b.store(false, std::memory_order_relaxed);
  }
  void concat_hashed_qgram_vectors(HashedQgramVector<sizeof_unit>
                                     *hashed_qgrams_vector) noexcept
  {
    size_t current_idx = 0;
    size_t max_size_idx = 0;
    size_t max_size = 0;
    size_t total_number_of_hashed_qgrams = 0;
    for (auto const &mv : table)
    {
      total_number_of_hashed_qgrams += mv.size();
      if (max_size < mv.size())
      {
        max_size = mv.size();
        max_size_idx = current_idx;
      }
      current_idx++;
    }
    if (table.size() > 1)
    {
      table[max_size_idx].reserve(total_number_of_hashed_qgrams);
      for (size_t idx = 0; idx < table.size(); idx++)
      {
        if (idx != max_size_idx)
        {
          table[max_size_idx].insert(table[max_size_idx].end(),
                                     table[idx].begin(),
                                     table[idx].end());
        }
      }
    }
    *hashed_qgrams_vector = std::move(table[max_size_idx]);
  }
  [[nodiscard]] size_t count_all_qgrams_get(void) const noexcept
  {
    size_t total_count_all_qgrams = 0;
    for (auto count : count_all_qgrams)
    {
      total_count_all_qgrams += count;
    }
    return total_count_all_qgrams;
  }
  [[nodiscard]] bool has_wildcards_get(void) const noexcept
  {
    size_t total_has_wildcards = false;
    for (const auto &hw : has_wildcards)
    {
      total_has_wildcards = total_has_wildcards || hw;
    }
    return total_has_wildcards;
  }
  ~HashedQgramVectorTable(void) = default;
};

template<int sizeof_unit,class HashIterator>
static void append_hashed_qgrams_threaded(size_t thread_id,
                                          size_t task_num,
                                          const GttlMultiseq &multiseq,
                                          size_t qgram_length,
                                          size_t window_size,
                                          uint64_t hash_mask,
                                          const GttlBitPacker<sizeof_unit,3>
                                            &hashed_qgram_packer,
                                          HashedQgramVectorTable<sizeof_unit>
                                            *hashed_qgram_vector_table)
{
  size_t this_count;
  size_t this_has_wildcards;
  std::tie(this_count,this_has_wildcards)
    = append_minimizers<sizeof_unit,HashIterator>
                       (qgram_length,
                        window_size,
                        hash_mask,
                        hashed_qgram_packer,
                        &hashed_qgram_vector_table->table[thread_id],
                        multiseq.sequence_ptr_get(task_num),
                        multiseq.sequence_length_get(task_num),
                        task_num);
  hashed_qgram_vector_table->count_all_qgrams[thread_id] += this_count;
  if(this_has_wildcards)
  {
    hashed_qgram_vector_table->has_wildcards[thread_id]
      .store(true, std::memory_order_relaxed);
  }
}

template<int sizeof_unit,class HashIterator>
class HashedQgramsGeneric
{
  struct Iterator
  {
    struct DecodedHashedQgram
    {
      uint64_t hash_value;
      size_t sequence_number,
             startpos;
      DecodedHashedQgram(uint64_t _hash_value,
                     size_t _sequence_number,
                     size_t _startpos)
       : hash_value(_hash_value)
       , sequence_number(_sequence_number)
       , startpos(_startpos)
      {}
      bool operator == (const DecodedHashedQgram &other) const noexcept
      {
        return hash_value == other.hash_value;
      }
      bool operator < (const DecodedHashedQgram &other) const noexcept
      {
        return hash_value < other.hash_value;
      }
    };
    const HashedQgramVector<sizeof_unit> &hashed_qgram_vector;
    const GttlBitPacker<sizeof_unit,3> &hashed_qgram_packer;
    size_t current_idx;
    Iterator(const HashedQgramVector<sizeof_unit> &_hashed_qgram_vector,
             const GttlBitPacker<sizeof_unit,3> &_hashed_qgram_packer,
             size_t _current_idx)
     : hashed_qgram_vector(_hashed_qgram_vector)
     , hashed_qgram_packer(_hashed_qgram_packer)
     , current_idx(_current_idx)
    {}
    public:
    const DecodedHashedQgram operator*(void) const noexcept
    {
      assert(current_idx < hashed_qgram_vector.size());
      return DecodedHashedQgram(hashed_qgram_vector[current_idx].
                                  template decode_at<0>(hashed_qgram_packer),
                                static_cast<size_t>(
                                  hashed_qgram_vector[current_idx].
                                    template decode_at<1>(hashed_qgram_packer)),
                                static_cast<size_t>(
                                  hashed_qgram_vector[current_idx].
                                    template decode_at<2>(hashed_qgram_packer))
                               );
    }
    auto& operator++() /* prefix increment*/
    {
      current_idx++;
      return *this;
    }
    bool operator != (const Iterator& other) const noexcept
    {
      return current_idx != other.current_idx;
    }
  };
  public:
  static constexpr int sizeof_unit_const = sizeof_unit;
  static constexpr bool handle_both_strands = HashIterator::handle_both_strands;
  const GttlMultiseq &multiseq;
  private:
  HashedQgramVector<sizeof_unit> hashed_qgram_vector;
  bool has_wildcards;
  int hashbits;
  size_t count_all_qgrams;
  size_t qgram_length;
  GttlBitPacker<sizeof_unit,3> hashed_qgram_packer;
  [[nodiscard]]
  size_t match_pos_pair_determine_length(uint64_t first_code,
                                         size_t start_range) const noexcept
  {
    size_t idx;

    for (idx = start_range + 1;
         idx < hashed_qgram_vector.size() &&
         first_code == this->hash_value_get(idx);
         idx++)
      /* Nothing */;
    return idx - start_range;
  }
  size_t remove_replicates_inplace(size_t max_replicates)
  {
    assert(max_replicates > 0);
    if (hashed_qgram_vector.size() <= 1)
    {
      return 0;
    }
    size_t r_idx = 0;
    size_t w_idx = 0;
    while (r_idx < hashed_qgram_vector.size())
    {
      const uint64_t first_code = this->hash_value_get(r_idx);
      const size_t replicate_len
        = match_pos_pair_determine_length(first_code,r_idx);
      assert(w_idx <= r_idx);
      if (replicate_len <= max_replicates)
      {
        if (w_idx < r_idx)
        {
          for (size_t idx = 0; idx < replicate_len; idx++)
          {
            hashed_qgram_vector[w_idx+idx] = hashed_qgram_vector[r_idx+idx];
          }
        }
        w_idx += replicate_len;
      }
      r_idx += replicate_len;
    }
    hashed_qgram_vector.resize(w_idx);
    return r_idx - w_idx;
  }
  public:
  static constexpr const bool possible_false_positive_matches
    = HashIterator::possible_false_positive_matches;
  [[nodiscard]] size_t size(void) const noexcept
  {
    return hashed_qgram_vector.size();
  }
  HashedQgramsGeneric(const GttlMultiseq &_multiseq,
                      size_t number_of_threads,
                      size_t _qgram_length,
                      size_t window_size,
                      int _hashbits,
                      bool sort_by_hashvalue,
                      bool at_constant_distance,
                      size_t max_replicates,
                      std::vector<std::string> *log_vector)
    : multiseq(_multiseq)
    , has_wildcards(false)
    , hashbits(_hashbits)
    , count_all_qgrams(0)
    , qgram_length(_qgram_length)
    , hashed_qgram_packer(GttlBitPacker<sizeof_unit,3>(
                            {_hashbits,
                             multiseq.sequences_number_bits_get(),
                             multiseq.sequences_length_bits_get()}))
  {
    assert(hashbits != -1);
    RunTimeClass rt_collect{};
    assert(number_of_threads >= 1);
    if (log_vector != nullptr)
    {
      log_vector->push_back(std::string("kmer_size\t") +
                            std::to_string(qgram_length));
      log_vector->push_back(std::string("window_size\t") +
                            std::to_string(window_size));
    }
    const uint64_t hash_mask = gttl_bits2maxvalue<uint64_t>(hashbits);
    if (number_of_threads == 1)
    {
      for (size_t seqnum = 0; seqnum < multiseq.sequences_number_get();
           seqnum++)
      {
        size_t this_count;
        size_t this_has_wildcards;
        std::tie(this_count,this_has_wildcards)
          = (at_constant_distance
               ? append_constant_distance_hashed_qgrams
                                  <sizeof_unit,HashIterator>
               : append_minimizers<sizeof_unit,HashIterator>)
                                  (qgram_length,
                                   window_size,
                                   hash_mask,
                                   hashed_qgram_packer,
                                   &hashed_qgram_vector,
                                   multiseq.sequence_ptr_get(seqnum),
                                   multiseq.sequence_length_get(seqnum),
                                   seqnum);
        count_all_qgrams += this_count;
        has_wildcards = has_wildcards || this_has_wildcards;
      }
    } else
    {
      assert(!at_constant_distance);
      HashedQgramVectorTable<sizeof_unit>
        hashed_qgram_vector_table(number_of_threads);
      GttlThreadPoolVar(number_of_threads,
                        multiseq.sequences_number_get(),
                        append_hashed_qgrams_threaded<sizeof_unit,HashIterator>,
                        multiseq,
                        qgram_length,
                        window_size,
                        hash_mask,
                        hashed_qgram_packer,
                        &hashed_qgram_vector_table);
      RunTimeClass rt_concat{};
      hashed_qgram_vector_table
        .concat_hashed_qgram_vectors(&hashed_qgram_vector);
      count_all_qgrams = hashed_qgram_vector_table.count_all_qgrams_get();
      has_wildcards = hashed_qgram_vector_table.has_wildcards_get();
      if (log_vector != nullptr)
      {
        log_vector->push_back(rt_concat.
                              to_string("concatenate hashed kmers vectors"));
      }
    }
    if (log_vector != nullptr)
    {
      log_vector->push_back(rt_collect.to_string("collect hashed kmers"));
    }

    if (sort_by_hashvalue)
    {
      RunTimeClass rt_postprocess{};
      if constexpr (sizeof_unit == 8)
      {
        uint64_t * const ptr =
          reinterpret_cast<uint64_t *>(hashed_qgram_vector.data());
        const Buckets<size_t> *const buckets
          = ska_lsb_radix_sort<size_t>(hashbits,ptr,hashed_qgram_vector.size());
        delete buckets;
      } else
      {
        const bool reversed_byte_order = not is_big_endian();
        ska_large_lsb_small_radix_sort(sizeof_unit,
                                       hashbits,
                                       reinterpret_cast<uint8_t *>
                                         (hashed_qgram_vector.data()),
                                       hashed_qgram_vector.size(),
                                       reversed_byte_order);
      }
      if (log_vector != nullptr)
      {
        log_vector->push_back(rt_postprocess.
                              to_string("sort hashed kmers by hash value"));
      }
      if (max_replicates > 0)
      {
        rt_postprocess.reset();
        const size_t orig_size = hashed_qgram_vector.size();
        const size_t removed = remove_replicates_inplace(max_replicates);
        if (log_vector != nullptr)
        {
          const double percentage
            = 100 * static_cast<double>(removed)/orig_size;
          auto t_msg = std::format("removed {} replicates with more than {} "
                                   "occurrences ({:.2f}% of all {} minimizers)",
                                   removed,
                                   max_replicates,
                                   percentage,
                                   orig_size);
          log_vector->push_back(rt_postprocess.to_string(t_msg));
        }
      }
    }
    if (log_vector != nullptr)
    {
      log_vector->push_back(std::string("number of hashed kmers\t") +
                            std::to_string(this->size()));
      const double density
        = static_cast<double>(size())/this->count_all_qgrams_get();
      log_vector->push_back(std::format("hashed kmers density\t{:.2f}",
                                        density));
      const double space_in_mega_bytes = mega_bytes(this->size() * sizeof_unit);
      log_vector->push_back(std::format("SPACE\thashed kmers (MB)\t{:.0f}",
                                        std::ceil(space_in_mega_bytes)));
    }
  }
  /* Construktor to be instantiated by Vector of bitpacker (Cached IO) */
  HashedQgramsGeneric(const GttlMultiseq &_multiseq,
                      size_t _qgram_length,
                      int _hashbits,
                      size_t _count_all_qgrams,
                      bool _has_wildcards,
                      const GttlBitPacker<sizeof_unit,3>
                        &_hashed_qgram_packer,
                      HashedQgramVector<sizeof_unit> &&_hashed_qgram_vector)
    : multiseq(_multiseq)
    , hashed_qgram_vector(std::move(_hashed_qgram_vector))
    , has_wildcards(_has_wildcards)
    , hashbits(_hashbits)
    , count_all_qgrams(_count_all_qgrams)
    , qgram_length(_qgram_length)
    , hashed_qgram_packer(_hashed_qgram_packer)
  { }
  [[nodiscard]] size_t count_all_qgrams_get(void) const noexcept
  {
    return count_all_qgrams;
  }
  [[nodiscard]] bool sequence_has_wildcards(void) const noexcept
  {
    return has_wildcards;
  }
  [[nodiscard]] uint64_t hash_value_get(size_t idx) const noexcept
  {
    assert(idx < size());
    return hashed_qgram_vector[idx].template decode_at<0>(hashed_qgram_packer);
  }
  [[nodiscard]] std::map<size_t, size_t> hash_value_run_statistics(void) const
  {
    std::map<size_t, size_t> count_runs;
    if (hashed_qgram_vector.size() > 0)
    {
      uint64_t previous_hash_value = hash_value_get(0);
      size_t run_length = size_t(1);
      for (size_t idx = 1; idx < hashed_qgram_vector.size(); idx++)
      {
        const uint64_t this_hash_value = hash_value_get(idx);
        if (previous_hash_value != this_hash_value)
        {
          count_runs[run_length]++;
          run_length = size_t(1);
          previous_hash_value = this_hash_value;
        } else
        {
          run_length++;
        }
      }
      count_runs[run_length]++;
      size_t total_runs = 0;
      for (auto &[r, c] : count_runs)
      {
        total_runs += r * c;
      }
      if (total_runs != hashed_qgram_vector.size())
      {
        throw std::runtime_error(
          std::format("sum of runs = {} != {} = number of minimizers",
                      total_runs, hashed_qgram_vector.size()));
      }
    }
    return count_runs;
  }
  [[nodiscard]] size_t sequence_number_get(size_t idx) const noexcept
  {
    assert(idx < size());
    return static_cast<size_t>(hashed_qgram_vector[idx]
                               .template decode_at<1>(hashed_qgram_packer));
  }
  [[nodiscard]] size_t startpos_get(size_t idx) const noexcept
  {
    assert(idx < size());
    return static_cast<size_t>(hashed_qgram_vector[idx]
                               .template decode_at<2>(hashed_qgram_packer));
  }
  [[nodiscard]] int packer_bit_group_size_get(int idx) const noexcept
  {
    return hashed_qgram_packer.bit_group_size_get(idx);
  }

  [[nodiscard]] size_t qgram_length_get(void) const noexcept
  {
    return qgram_length;
  }
  /* the following method has a side effect */
  void show(size_t offset = 0) const noexcept
  {
    printf("# Hash\tSeqNr\tStart\n");
    for (size_t idx = 0; idx < size(); idx++)
    {
      const uint64_t hash_value = hash_value_get(idx);
      const size_t sequence_number = sequence_number_get(idx);
      const size_t startpos = startpos_get(idx);
      printf("%zu\t%zu\t%zu\n",static_cast<size_t>(hash_value),
             sequence_number + offset,startpos);
    }
  }
  [[nodiscard]] Iterator begin(void) const noexcept
  {
    return Iterator (hashed_qgram_vector,hashed_qgram_packer,0);
  }
  [[nodiscard]] Iterator end(void) const noexcept
  {
    return Iterator (hashed_qgram_vector,hashed_qgram_packer,
                     hashed_qgram_vector.size());
  }
  /* save index to file */
  [[nodiscard]] bool save_cache(const std::string &cache_path,
                                uint64_t ref_window_size,
                                bool at_constant_distance_flag) const
  {
    std::ofstream os(cache_path,std::ios::binary | std::ios::trunc);
    if (!os)
    {
      return false;
    }
    HashedQgramsCacheHeader header{};
    std::memcpy(header.magic,hashed_qgrams_cache_magic,
                sizeof hashed_qgrams_cache_magic);
    header.version = hashed_qgrams_cache_version;
    header.sizeof_unit = static_cast<uint8_t>(sizeof_unit);
    header.handle_both_strands
      = HashIterator::handle_both_strands ? 1 : 0;
    header.possible_false_positive
      = HashIterator::possible_false_positive_matches ? 1 : 0;
    header.at_constant_distance = at_constant_distance_flag ? 1 : 0;
    header.has_wildcards = has_wildcards ? 1 : 0;
    header.qgram_length = static_cast<uint64_t>(qgram_length);
    header.hashbits = static_cast<int64_t>(hashbits);
    header.count_all_qgrams = static_cast<uint64_t>(count_all_qgrams);
    header.ref_window_size = ref_window_size;
    header.sequences_number_bits
      = static_cast<uint64_t>(multiseq.sequences_number_bits_get());
    header.sequences_length_bits
      = static_cast<uint64_t>(multiseq.sequences_length_bits_get());
    header.hashed_qram_vector_size
      = static_cast<uint64_t>(hashed_qgram_vector.size());
    os.write(reinterpret_cast<const char *>(&header),sizeof header);
    if (!os)
    {
      return false;
    }
    if (!hashed_qgram_vector.empty())
    {
      os.write(reinterpret_cast<const char *>(hashed_qgram_vector.data()),
               hashed_qgram_vector.size()
                 * sizeof (BytesUnit<sizeof_unit,3>));
      if (!os)
      {
        return false;
      }
    }
    return true;
  }

  /* load index from cache wieder, fails if parameters are incorrect */
  static HashedQgramsGeneric /* read header, vector, bitpacker */
    load_cache(const std::string &cache_path,
               const GttlMultiseq &multiseq,
               size_t expected_qgram_length,
               int expected_hashbits,
               uint64_t ref_window_size,
               bool at_constant_distance_flag)
  {
    std::ifstream is(cache_path,std::ios::binary);
    if (!is)
    {
      throw std::runtime_error(std::format("cache file {}: cannot open",
                                           cache_path));
    }
    HashedQgramsCacheHeader header{};
    is.read(reinterpret_cast<char *>(&header),sizeof header);
    if (!is)
    {
      throw std::runtime_error(std::format("cache file {}: cannot read header",
                                           cache_path));
    }
    if (std::memcmp(header.magic,hashed_qgrams_cache_magic,
                    sizeof hashed_qgrams_cache_magic) != 0)
    {
      throw std::runtime_error(std::format("cache file {}: invalid magic "
                                           "header",cache_path));
    }
    if (header.version != hashed_qgrams_cache_version)
    {
      throw std::runtime_error(std::format("cache file {}: unsupported cache "
                                           "version {}",
                                           cache_path,header.version));
    }
    if (header.sizeof_unit != static_cast<uint8_t>(sizeof_unit))
    {
      throw std::runtime_error(std::format("cache file {}: incompatible "
                                           "sizeof_unit ({} != {})",
                                           cache_path,
                                           sizeof_unit,
                                           header.sizeof_unit));
    }
    const uint8_t expected_handle_both_strands
      = HashIterator::handle_both_strands ? 1 : 0;
    if (header.handle_both_strands != expected_handle_both_strands)
    {
      throw std::runtime_error(std::format("cache file {}: reverse-complement "
                                           " handling differs from cached "
                                           "version",cache_path));
    }
    const uint8_t expected_possible_false_positive
      = HashIterator::possible_false_positive_matches ? 1 : 0;
    if (header.possible_false_positive != expected_possible_false_positive)
    {
      throw std::runtime_error(std::format("cache file {}: hash method "
                                           "not compatible with cached version",
                                           cache_path));
    }
    if (header.at_constant_distance
          != static_cast<uint8_t>(at_constant_distance_flag ? 1 : 0))
    {
      throw std::runtime_error(std::format("cache file {}: constant-distance "
                                           "flag not compatible with cached "
                                           "version", cache_path));
    }
    if (header.qgram_length != expected_qgram_length)
    {
      throw std::runtime_error(std::format("cache file {}: kmer size {} is "
                                           "not comptable with cached version "
                                           "which was created for k={}",
                                           cache_path,
                                           expected_qgram_length,
                                           header.qgram_length));
    }
    if (header.hashbits != expected_hashbits)
    {
      throw std::runtime_error(std::format("cache file {}: hash_bits={} is "
                                           "not compatible with {} as used in "
                                           "the cached version",
                                           cache_path,
                                           expected_hashbits,
                                           header.hashbits));
    }
    if (header.ref_window_size != ref_window_size)
    {
      throw std::runtime_error(std::format("cache file {}: "
                                           "window_size={} is not compatible "
                                           "with {} as used in the cached "
                                           "version",
                                           cache_path,
                                           ref_window_size,
                                           header.ref_window_size));
    }
    if (std::cmp_not_equal(header.sequences_number_bits,
                           multiseq.sequences_number_bits_get()))
    {
      throw std::runtime_error(std::format("cache file {}: sequence number "
                                           "bits {} of reference sequence is "
                                           "not compatible with {} as used in "
                                           "the cached verstion",
                                           cache_path,
                                           multiseq.sequences_number_bits_get(),
                                           header.sequences_number_bits));
    }
    if (std::cmp_not_equal(header.sequences_length_bits,
                           multiseq.sequences_length_bits_get()))
    {
      throw std::runtime_error(std::format("cache file {}: sequence length "
                                           "bits {} of reference sequence is "
                                           "not compatible with {} as used in "
                                           "the cached verstion",
                                           cache_path,
                                           multiseq.sequences_length_bits_get(),
                                           header.sequences_length_bits));
    }
    if (header.hashed_qram_vector_size
          > static_cast<uint64_t>(SIZE_MAX / sizeof (BytesUnit<sizeof_unit,3>)))
    {
      throw std::runtime_error(std::format("cache file {}: hashed qgram vector "
                                           "is too large",
                                           cache_path));
    }
    HashedQgramVector<sizeof_unit>
      hashed_qgrams(static_cast<size_t>(header.hashed_qram_vector_size));
    if (!hashed_qgrams.empty())
    {
      is.read(reinterpret_cast<char *>(hashed_qgrams.data()),
              hashed_qgrams.size() * sizeof (BytesUnit<sizeof_unit,3>));
      if (!is)
      {
        throw std::runtime_error(std::format("cache file {}: truncated hashed "
                                             "qgram data",
                                             cache_path));
      }
    }
    const GttlBitPacker<sizeof_unit,3> packer(
      {expected_hashbits,
       multiseq.sequences_number_bits_get(),
       multiseq.sequences_length_bits_get()});
    return HashedQgramsGeneric(multiseq,
                               expected_qgram_length,
                               expected_hashbits,
                               static_cast<size_t>(header.count_all_qgrams),
                               header.has_wildcards != 0,
                               packer,
                               std::move(hashed_qgrams));
  }
};
#endif
