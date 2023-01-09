#ifndef HASHED_QGRAMS_HPP
#define HASHED_QGRAMS_HPP
#include <cstddef>
#include <cstdint>
#include <cassert>
#include <cstdlib>
#include <cstdbool>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <ios>
#include <deque>
#include <type_traits>

#include "utilities/ska_lsb_radix_sort.hpp"
#include "utilities/mathsupport.hpp"
#include "utilities/runtime_class.hpp"
#include "utilities/bytes_unit.hpp"
#include "utilities/is_big_endian.hpp"
#include "utilities/str_format.hpp"
#include "threading/thread_pool_var.hpp"
#include "sequences/gttl_seq_iterator.hpp"
#include "sequences/char_range.hpp"
#include "sequences/char_finder.hpp"
#include "sequences/gttl_multiseq.hpp"

template<int sizeof_unit>
using HashedQgramVector = std::vector<BytesUnit<sizeof_unit,3>>;

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
      StrFormat msg("cannot find minimizer %lu",static_cast<size_t>(min_hash));
      throw msg.str();
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
  NucleotideRanger ranger(sequence,seqlen);
  HashedQgramVector<sizeof_unit> palindromic_vector{};
  for (auto const &&range : ranger)
  {
    const size_t this_length = std::get<1>(range);
    has_wildcards = has_wildcards || this_length < seqlen;
    if (this_length < minseqlen_to_process)
    {
      continue;
    }
    size_t seqpos = std::get<0>(range);
    const char *seqptr = sequence + seqpos;
    HashIterator qgiter(qgram_length, seqptr, this_length);
    count_all_qgrams += (this_length - qgram_length - 1);

#ifdef VALIDATE_MINIMIZER
    HashedQgramVector<sizeof_unit> all_hashed_qgrams{};
    size_t search_start = 0;
#endif
    bool front_was_moved = false;
    for (auto const &&code_pair : qgiter)
    {
      uint64_t this_hash = std::get<0>(code_pair) & hash_mask;
      size_t stored_seqnum, stored_seqpos;
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
            BytesUnit<sizeof_unit,3>
              palindromic_hashed_qgram_rc(hashed_qgram_packer,
                                       {this_hash,
                                        static_cast<uint64_t>(seqnum+1),
                                        static_cast<uint64_t>(seqlen -
                                                              (seqpos +
                                                               qgram_length))});
            palindromic_vector.emplace_back(palindromic_hashed_qgram_rc);
            BytesUnit<sizeof_unit,3>
              palindromic_hashed_qgram_fwd(hashed_qgram_packer,
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
      BytesUnit<sizeof_unit,3>
          current_hashed_qgram(hashed_qgram_packer,
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
        uint64_t min_pos_in_window
          = window_deque.front().template decode_at<2>(hashed_qgram_packer);
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
    std::cout << "# add " << palindromic_vector.size() << " " << qgram_length
              << "-mer DNA-palindromes to minimizer_vector" << std::endl;
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
  NucleotideRanger ranger(sequence,seqlen);
  for (auto const &&range : ranger)
  {
    const size_t this_length = std::get<1>(range);
    has_wildcards = has_wildcards || this_length < seqlen;
    if (this_length < minseqlen_to_process)
    {
      continue;
    }
    size_t seqpos = std::get<0>(range);
    const char *seqptr = sequence + seqpos;
    HashIterator qgiter(qgram_length, seqptr, this_length);
    count_all_qgrams += (this_length - qgram_length - 1);
    size_t steps = 0;
    for (auto const &&code_pair : qgiter)
    {
      if (steps == 0)
      {
        uint64_t this_hash = std::get<0>(code_pair) & hash_mask;
        size_t stored_seqnum, stored_seqpos;
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
        BytesUnit<sizeof_unit,3>
          current_hashed_qgram(hashed_qgram_packer,
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
  std::vector<bool> has_wildcards;
  std::vector<HashedQgramVector<sizeof_unit>> table;
  HashedQgramVectorTable(size_t number_of_threads)
    : count_all_qgrams(number_of_threads,0)
    , has_wildcards(number_of_threads,false)
    , table(number_of_threads,HashedQgramVector<sizeof_unit>{})
  {
  }
  void concat_hashed_qgram_vectors(HashedQgramVector<sizeof_unit>
                                     *hashed_qgrams_vector) noexcept
  {
    size_t current_idx = 0, max_size_idx = 0, max_size = 0,
           total_number_of_hashed_qgrams = 0;
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
  size_t count_all_qgrams_get(void) const noexcept
  {
    size_t total_count_all_qgrams = 0;
    for (auto count : count_all_qgrams)
    {
      total_count_all_qgrams += count;
    }
    return total_count_all_qgrams;
  }
  bool has_wildcards_get(void) const noexcept
  {
    size_t total_has_wildcards = false;
    for (auto hw : has_wildcards)
    {
      total_has_wildcards = total_has_wildcards || hw;
    }
    return total_has_wildcards;
  }
  ~HashedQgramVectorTable(void) {}
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
  size_t this_count, this_has_wildcards;
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
  hashed_qgram_vector_table->has_wildcards[thread_id] =
    hashed_qgram_vector_table->has_wildcards[thread_id] || this_has_wildcards;
}

template<int sizeof_unit,class HashIterator>
class HashedQgramsGeneric
{
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
  public:
  static constexpr const bool possible_false_positive_matches
    = HashIterator::possible_false_positive_matches;
  size_t size(void) const noexcept
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
                      std::vector<std::string> *log_vector)
    : multiseq(_multiseq)
    , hashed_qgram_vector({})
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
        size_t this_count, this_has_wildcards;
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
      log_vector->push_back(std::string("number of hashed kmers\t") +
                            std::to_string(this->size()));
      const double density
        = static_cast<double>(size())/this->count_all_qgrams_get();
      StrFormat s_density("hashed kmers density\t%.2f",density);
      log_vector->push_back(s_density.str());
      StrFormat s_space("SPACE\thashed kmers (MB)\t%lu",
                        static_cast<unsigned long>(std::ceil(
                                                   mega_bytes(this->size()
                                                              * sizeof_unit))));
      log_vector->push_back(s_space.str());
      log_vector->push_back(rt_collect.to_string("collect hashed kmers"));
    }
    if (sort_by_hashvalue)
    {
      RunTimeClass rt_sort{};
      if constexpr (sizeof_unit == 8)
      {
        Buckets<size_t> *buckets
          = ska_lsb_radix_sort<size_t>(hashbits,
                                       reinterpret_cast<uint64_t *>
                                         (hashed_qgram_vector.data()),
                                       hashed_qgram_vector.size());
        delete buckets;
      } else
      {
        const bool reversed_byte_order = is_big_endian() ? false : true;
        ska_large_lsb_small_radix_sort(sizeof_unit,
                                       hashbits,
                                       reinterpret_cast<uint8_t *>
                                         (hashed_qgram_vector.data()),
                                       hashed_qgram_vector.size(),
                                       reversed_byte_order);
      }
      if (log_vector != nullptr)
      {
        log_vector->push_back(rt_sort.
                              to_string("sort hashed kmers by hash value"));
      }
    }
  }
  size_t count_all_qgrams_get(void) const noexcept
  {
    return count_all_qgrams;
  }
  bool sequence_has_wildcards(void) const noexcept
  {
    return has_wildcards;
  }
  uint64_t hash_value_get(size_t idx) const noexcept
  {
    assert(idx < size());
    return hashed_qgram_vector[idx].template decode_at<0>(hashed_qgram_packer);
  }
  size_t sequence_number_get(size_t idx) const noexcept
  {
    assert(idx < size());
    return static_cast<size_t>(hashed_qgram_vector[idx]
                               .template decode_at<1>(hashed_qgram_packer));
  }
  size_t startpos_get(size_t idx) const noexcept
  {
    assert(idx < size());
    return static_cast<size_t>(hashed_qgram_vector[idx]
                               .template decode_at<2>(hashed_qgram_packer));
  }
  int packer_bit_group_size_get(int idx) const noexcept
  {
    return hashed_qgram_packer.bit_group_size_get(idx);
  }

  size_t qgram_length_get(void) const noexcept
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
      printf("%lu\t%lu\t%lu\n",static_cast<unsigned long>(hash_value),
             sequence_number + offset,startpos);
    }
  }
};
#endif
