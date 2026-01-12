#ifndef GTTL_MINIMIZER_GENERATOR_HPP
#define GTTL_MINIMIZER_GENERATOR_HPP

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <optional>
#include <string_view>
#include <vector>
#include "sequences/char_finder.hpp"
#include "sequences/char_range.hpp"

// MinimizerValueClass is a std::tuple<uint64_t, size_t, size_t>
// in our testsuite code. It may instead be a bitpacker or similar type,
// given that it supports construction from a 3-tuple
// (hash_value, seqnum, position)
// and can be copied, moved and stored in an STL-container.
template <class HashIterator, class MinimizerValueClass>
class GttlHashedMinimizerGenerator
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

  // We use std::optionals here so that these variables are NEVER
  // default-initialized. This is necessary because there is no default
  // constructor for them.
  std::optional<decltype(qgiter->begin())> current_it;
  std::optional<decltype(qgiter->end())> current_it_end;

  bool front_was_moved;
  bool using_palindromic;

  public:

  explicit GttlHashedMinimizerGenerator(size_t _qgram_length,
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
  GttlHashedMinimizerGenerator(const GttlHashedMinimizerGenerator&)
    = delete;
  GttlHashedMinimizerGenerator& operator=(const GttlHashedMinimizerGenerator&)
    = delete;
  GttlHashedMinimizerGenerator(GttlHashedMinimizerGenerator&&)
    = delete;
  GttlHashedMinimizerGenerator& operator=(GttlHashedMinimizerGenerator&&)
    = delete;

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
          // "emplace" in this case simply populates a std::optional
          current_it.emplace(qgiter->begin());
          current_it_end.emplace(qgiter->end());

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
      while (current_it != current_it_end)
      {
        const auto &code_pair = *(current_it.value());
        ++current_it.value();

        // NOLINTNEXTLINE(misc-const-correctness)
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
          // NOLINTNEXTLINE(misc-const-correctness)
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
    GttlHashedMinimizerGenerator* generator;
    bool is_end;
    public:
    explicit Iterator(GttlHashedMinimizerGenerator* _generator,
                      bool end = false)
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





// KeyBuilder must have an operator() that produces
// objects which can be compared using a < relation.
// An example might be our own UHSKeyBuilder or a functor
// class that simply calls std::hash.
template <class MinimizerValueClass, class KeyBuilder>
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
  const std::string_view sequence;
  const size_t seqlen;
  const size_t minseqlen_to_process;

  // Internal state
  static constexpr const char_finder::NucleotideFinder unw_nucleotide_finder{};

  std::deque<MinimizerValueClass> window_deque;

  size_t this_length;
  size_t seqpos;

  bool front_was_moved;
  size_t last_emitted_pos = SIZE_MAX;

  const KeyBuilder key_builder{};

  public:

  explicit GttlMinimizerGenerator(size_t _qgram_length,
                                  size_t _window_size,
                                  std::string_view _sequence,
                                  KeyBuilder kb,
                                  MinimizerValueClass* _out = nullptr,
                                  bool _is_end = false)
    : out(_out == nullptr ? &default_buffer : _out)
    , is_end(_is_end)
    , qgram_length(_qgram_length)
    , window_size(_window_size)
    , sequence(_sequence)
    , seqlen(_sequence.size())
    , minseqlen_to_process(window_size + qgram_length - 1)
    , this_length(0)
    , seqpos(0)
    , front_was_moved(false)
    , key_builder(kb)
  {
    if (seqlen < minseqlen_to_process)
    {
      is_end = true;
    }
  }

  // Delete copy/move constructor & assignment operator
  GttlMinimizerGenerator(const GttlMinimizerGenerator&) = delete;
  GttlMinimizerGenerator& operator=(const GttlMinimizerGenerator&) = delete;
  GttlMinimizerGenerator(GttlMinimizerGenerator&&) = delete;
  GttlMinimizerGenerator& operator=(GttlMinimizerGenerator&&) = delete;

  bool advance(void)
  {
    if (is_end) [[unlikely]]
    {
      return false;
    }

    while (seqpos + qgram_length <= seqlen)
    {
      const bool out_of_window =
        seqpos >= window_size and
        std::get<1>(window_deque.front()) <= seqpos - window_size;

      if (not window_deque.empty() and out_of_window)
      {
        window_deque.pop_front();
      }

      const auto key = key_builder(sequence.substr(seqpos, qgram_length));

      while (not window_deque.empty() and
             key < std::get<0>(window_deque.back()))
      {
        window_deque.pop_back();
      }

      window_deque.emplace_back(MinimizerValueClass({key,
                                                     seqpos}));

      if (seqpos >= window_size - 1) [[unlikely]]
      {
        const size_t current_min = std::get<1>(window_deque.front());
        if (last_emitted_pos != current_min)
        {
          *out = window_deque.front();
          last_emitted_pos = current_min;
          ++seqpos;
          return true;
        }
      }

      ++seqpos;
    }

    is_end = true;
    return false;
  }

  void reset(void)
  {
    seqpos = 0;
    window_deque.clear();
    is_end = seqlen < minseqlen_to_process;
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
