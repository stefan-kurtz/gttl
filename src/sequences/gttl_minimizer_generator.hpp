#ifndef GTTL_MINIMIZER_GENERATOR_HPP
#define GTTL_MINIMIZER_GENERATOR_HPP

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <utility>
#include <vector>
#include "sequences/char_finder.hpp"
#include "sequences/char_range.hpp"

template <class HashIterator, class MinimizerValueClass>
class GttlMinimizerGenerator
{
  private:
  // Default output buffer, in case the caller does not provide one
  MinimizerValueClass default_buffer;
  // output-pointer towards a minimizer value where the
  // current minimizer should be stored.
  MinimizerValueClass* out;
  // Whether all minimizers have been generated
  bool is_end;

  // Minimizer parameters (these are constant for the lifetime of the generator)
  const size_t qgram_length;
  const size_t window_size;
  const uint64_t hash_mask;
  const char* sequence;
  const size_t seqlen;
  const size_t seqnum;
  const size_t minseqlen_to_process;

  // Internal state
  static constexpr const char_finder::NucleotideFinder unw_nucleotide_finder{};

  std::deque<MinimizerValueClass> window_deque;

  GttlCharRange<char_finder::NucleotideFinder,
                unw_nucleotide_finder,
                true,
                false> ranger;

  decltype(ranger.begin()) range_it;
  decltype(ranger.end()) range_end;

  std::vector<MinimizerValueClass> palindromic_vector;
  size_t palindromic_pos;

  const char* seqptr;
  size_t this_length;
  size_t seqpos;

  HashIterator* qgiter;
  bool front_was_moved;
  bool using_palindromic;

  public:

  explicit GttlMinimizerGenerator(size_t _qgram_length,
                                  size_t _window_size,
                                  uint64_t _hash_mask,
                                  const char* _sequence,
                                  size_t _seqlen,
                                  size_t _seqnum,
                                  MinimizerValueClass* _out = nullptr,
                                  bool _is_end = false)
    : out(_out == nullptr ? &default_buffer : _out)
    , is_end(_is_end)
    , qgram_length(_qgram_length)
    , window_size(_window_size)
    , hash_mask(_hash_mask)
    , sequence(_sequence)
    , seqlen(_seqlen)
    , seqnum(_seqnum)
    , minseqlen_to_process(window_size + qgram_length - 1)
    , ranger(sequence, seqlen)
    , range_it(ranger.begin())
    , range_end(ranger.end())
    , palindromic_pos(0)
    , seqptr(nullptr)
    , this_length(0)
    , seqpos(0)
    , qgiter(nullptr)
    , front_was_moved(false)
    , using_palindromic(false)
  {
    if constexpr (HashIterator::handle_both_strands)
    {
      if (seqnum % 2 == 1)
      {
        is_end = true;
      }
    }
  }

  // Delete copy/move constructor & assignment operator
  GttlMinimizerGenerator(const GttlMinimizerGenerator&) = delete;
  GttlMinimizerGenerator& operator=(const GttlMinimizerGenerator&) = delete;
  GttlMinimizerGenerator(GttlMinimizerGenerator&&) = delete;
  GttlMinimizerGenerator& operator=(GttlMinimizerGenerator&&) = delete;

  bool advance(void)
  {
    if (is_end)
    {
      return false;
    }

    // Flush palindromic minimizers after main pass
    if (using_palindromic)
    {
      if (palindromic_pos < palindromic_vector.size())
      {
        *out = palindromic_vector[palindromic_pos++];
        return true;
      }
      is_end = true;
      return false;
    }

    while (true)
    {
      // Initialize new ranger segment
      if (qgiter == nullptr)
      {
        while (range_it != range_end)
        {
          const auto &range = *range_it;
          ++range_it;

          this_length = std::get<1>(range);
          if (this_length < minseqlen_to_process)
          {
            continue;
          }

          seqpos = std::get<0>(range);
          seqptr = sequence + seqpos;

          qgiter = new HashIterator(qgram_length, seqptr, this_length);

          window_deque.clear();
          front_was_moved = false;
          break;
        }

        if (qgiter == nullptr)
        {
          if constexpr (HashIterator::handle_both_strands)
          {
            using_palindromic = true;
            return advance();
          }

          is_end = true;
          return false;
        }
      }

      // iterate qgrams
      for (auto it = qgiter->begin(), it_end = qgiter->end();
           it != it_end; ++it)
      {
        const auto &code_pair = *it;
        ++it;

        uint64_t this_hash = std::get<0>(code_pair) & hash_mask;

        size_t stored_seqnum;
        size_t stored_seqpos;

        if constexpr (HashIterator::handle_both_strands)
        {
          const uint64_t rc_hash = std::get<1>(code_pair) & hash_mask;

          if (rc_hash < this_hash)
          {
            this_hash = rc_hash;
            stored_seqnum = seqnum + 1;
            stored_seqpos = seqlen - seqpos - qgram_length;
          }else
          {
            if (rc_hash == this_hash)
            {
              palindromic_vector.emplace_back(
                MinimizerValueClass({this_hash,
                                     static_cast<uint64_t>(seqnum+1),
                                     static_cast<uint64_t>(seqlen
                                                           - seqpos
                                                           - qgram_length)}));
              palindromic_vector.emplace_back(
                MinimizerValueClass({this_hash,
                                     static_cast<uint64_t>(seqnum),
                                     static_cast<uint64_t>(seqpos)}));
            }
            stored_seqnum = seqnum;
            stored_seqpos = seqpos;
          }
        }
        else
        {
          stored_seqnum = seqnum;
          stored_seqpos = seqpos;
        }

        while (not window_deque.empty()
               and std::get<0>(window_deque.back()) > this_hash)
        {
          window_deque.pop_back();
        }

        window_deque.emplace_back(
          MinimizerValueClass({this_hash,
                               static_cast<uint64_t>(stored_seqnum),
                               static_cast<uint64_t>(stored_seqpos)}));

        if (seqpos >= window_size)
        {
          uint64_t min_pos = std::get<2>(window_deque.front());

          if constexpr (HashIterator::handle_both_strands)
          {
            if (std::get<1>(window_deque.front()) == seqnum + 1)
            {
              min_pos = seqlen - min_pos - qgram_length;
            }
          }

          if (min_pos <= seqpos - window_size)
          {
            window_deque.pop_front();
            front_was_moved = false;
          }

          if (not front_was_moved)
          {
            *out = window_deque.front();
            front_was_moved = true;
            ++seqpos;
            return true;
          }
        }else if(seqpos == window_size - 1)
        {
          *out = window_deque.front();
          front_was_moved = true;
          ++seqpos;
          return true;
        }
        ++seqpos;
      }
      delete qgiter;
      qgiter = nullptr;
      window_deque.clear();
    }
  }

  void reset(void)
  {
    ranger.reset();
    range_it = ranger.begin();
    palindromic_vector.clear();
    palindromic_pos = 0;
    qgiter = nullptr;
    is_end = false;
    using_palindromic = false;
  }

  class Iterator
  {
    private:
    GttlMinimizerGenerator* generator;
    bool is_end;
    public:
    explicit Iterator(GttlMinimizerGenerator* _generator, bool end = false)
      : generator(_generator)
      , is_end(end)
    {
      if (not end)
      {
        ++(*this);
      }
    }

    const MinimizerValueClass* operator * () const
    {
      return generator->out;
    }

    const Iterator& operator ++ (void)
    {
      assert(not is_end);
      if (not generator->advance())
      {
        is_end = true;
      }
      return *this;
    }

    bool operator == (const Iterator& other) const
    {
      return is_end == other.is_end and generator == other.generator;
    }

    bool operator != (const Iterator& other) const
    {
      return not (*this == other);
    }

  };

  Iterator begin(void)
  {
    return Iterator(this, false);
  }

  Iterator end(void)
  {
    return Iterator(this, true);
  }
};

#endif // GTTL_MINIMIZER_GENERATOR_HPP
