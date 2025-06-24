#ifndef ENUMERATE_MINIMIZER_HPP
#define ENUMERATE_MINIMIZER_HPP
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>
#include <deque>
#include "sequences/char_range.hpp"
#include "sequences/char_finder.hpp"

template<class HashIterator,class MinimizerProcessor,class MinimizerValueClass>
void inline enumerate_minimizer(size_t qgram_length,
                                size_t window_size,
                                uint64_t hash_mask,
                                MinimizerProcessor *minimizer_processor,
                                const char *sequence,
                                size_t seqlen,
                                size_t seqnum)
{
  static constexpr const char_finder::NucleotideFinder unw_nucleotide_finder{};
  if constexpr (HashIterator::handle_both_strands)
  {
    if (seqnum % 2 == 1) /* for processing the reverse complement
                            we skip every second sequence, as the
                            minimizers come from the original sequence */
    {
      return;
    }
  }
  const size_t minseqlen_to_process = window_size + qgram_length - 1;
  std::deque<MinimizerValueClass> window_deque{};
  GttlCharRange<char_finder::NucleotideFinder,unw_nucleotide_finder,
                true, false> ranger(sequence,seqlen);
  std::vector<MinimizerValueClass> palindromic_vector{};
  for (auto const &&range : ranger)
  {
    const size_t this_length = std::get<1>(range);
    if (this_length < minseqlen_to_process)
    {
      continue;
    }
    size_t seqpos = std::get<0>(range);
    const char *seqptr = sequence + seqpos;
    HashIterator qgiter(qgram_length, seqptr, this_length);

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
            MinimizerValueClass
              palindromic_hashed_qgram_rc({this_hash,
                                           static_cast<uint64_t>(seqnum+1),
                                           static_cast<uint64_t>(seqlen -
                                                                 (seqpos +
                                                                  qgram_length))
                                          });
            palindromic_vector.emplace_back(palindromic_hashed_qgram_rc);
            MinimizerValueClass
              palindromic_hashed_qgram_fwd({this_hash,
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
        if (std::get<0>(window_deque.back()) <= this_hash)
        {
          break;
        }
        window_deque.pop_back();
      }
      MinimizerValueClass current_hashed_qgram({this_hash,
                                           static_cast<uint64_t>(stored_seqnum),
                                           static_cast<uint64_t>(stored_seqpos)}
                                          );
      window_deque.emplace_back(current_hashed_qgram);

      // After the minimizer of the first window was found
      // At this point we are only looking if
      // the current minimizer drops out of the window
      if (seqpos >= window_size)
      {
        uint64_t min_pos_in_window = std::get<2>(window_deque.front());
        if constexpr (HashIterator::handle_both_strands)
        {
          if (std::get<1>(window_deque.front()) == seqnum + 1)
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
        assert(not window_deque.empty()); /* because we added
                                             current_hashed_qgram */
        /* check if new minimizer was found (i.e. front element which was
           not moved yet */
        if (not front_was_moved)
        {
          minimizer_processor->apply(sequence,window_deque.front());
          front_was_moved = true; /* we moved front element and do not
                                     want to do it again */
        }
      } else
      {
        // add minimizer of first window
        if (seqpos == window_size - 1)
        {
          minimizer_processor->apply(sequence,window_deque.front());
          front_was_moved = true;
        }
      }
      seqpos++;
    }
    window_deque.clear();
  }
  if constexpr (HashIterator::handle_both_strands)
  {
    for (auto && pq : palindromic_vector)
    {
      minimizer_processor->apply(sequence,pq);
    }
  }
}
#endif
