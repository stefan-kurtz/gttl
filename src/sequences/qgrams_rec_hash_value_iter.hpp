/*
  Copyright (c) 2022 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2022 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/
#ifndef QGRAMS_REC_HASH_VALUE_ITER_HPP
#define QGRAMS_REC_HASH_VALUE_ITER_HPP
#include <cstddef>
#include <cstdint>
#include <cassert>
#include <cmath>
#include <iterator>
#include <tuple>
#include <utility>
#include "utilities/cyclic_buffer.hpp"
#include "sequences/alphabet.hpp"
#include "sequences/max_qgram_length.hpp"
#include "sequences/complement_uint8.hpp"

/* Implementation of iterator class follows concept described in
   https://davidgorski.ca/posts/stl-iterators/ */

template<const char *_char_spec,
         uint8_t _undefined_rank,
         class QgramTransformer>
class QgramRecHashValueIterator
{
  public:
  static constexpr const bool possible_false_positive_matches
    = QgramTransformer::possible_false_positive_matches;
  static constexpr const GttlAlphabet<_char_spec,_undefined_rank> alphabet{};
  static constexpr const bool handle_both_strands = true;
  /* As we use the expression uses 3 - base to compute the
     complement of a rank for a base, the following must be checked */
  static_assert(alphabet.char_to_rank('A') == 0 &&
                alphabet.char_to_rank('C') == 1 &&
                alphabet.char_to_rank('G') == 2 &&
                alphabet.char_to_rank('T') == 3);
  private:

  static constexpr const size_t alpha_size = alphabet.size();
  using SequenceBaseType = char;
  using CyclicBuffer_uint8 = CyclicBuffer<uint8_t,MAX_QGRAM_LENGTH>;

  struct Iterator
  {
    private:
      const QgramTransformer &qgram_transformer;
      CyclicBuffer_uint8 &current_window;
      const SequenceBaseType *next_char_ptr;
      bool last_qgram_was_processed;
      const SequenceBaseType *end_of_sequence;
      uint64_t hash_value, compl_hash_value;
      uint8_t wildcards_in_qgram;
    public:
      using iterator_category = std::forward_iterator_tag;
      using value_type = SequenceBaseType;
      using difference_type = size_t;
      using pointer = const SequenceBaseType*;
      using reference = uint64_t&;

      /* Constructor for begin() */
      Iterator(const QgramTransformer &_qgram_transformer,
               size_t _qgram_length,
               CyclicBuffer_uint8 &_current_window,
               pointer _sequence,
               uint64_t first_hash_value,
               uint64_t first_compl_hash_value)
        : qgram_transformer(_qgram_transformer)
        , current_window(_current_window)
        , next_char_ptr(_sequence + _qgram_length - 1)
        , last_qgram_was_processed(true)
        , end_of_sequence(nullptr)
        , hash_value(first_hash_value)
        , compl_hash_value(first_compl_hash_value)
      {
      }
      /* Constructor for end() */
      Iterator(const QgramTransformer &_qgram_transformer,
               CyclicBuffer_uint8 &_current_window,
               pointer _end_of_sequence)
        : qgram_transformer(_qgram_transformer)
        , current_window(_current_window)
        , end_of_sequence(_end_of_sequence)
      {
      }
      std::pair<uint64_t,uint64_t> operator*()
      {
        if (!last_qgram_was_processed)
        {
          const uint8_t new_rank = alphabet.char_to_rank(*next_char_ptr);
          const uint8_t old_rank = current_window.shift(new_rank);
          hash_value = qgram_transformer.next_hash_value_get(old_rank,
                                                             hash_value,
                                                             new_rank);
          compl_hash_value
            = qgram_transformer.next_compl_hash_value_get(
                 complement_uint8(old_rank),
                 compl_hash_value,
                 complement_uint8(new_rank));
        }
        return {hash_value,compl_hash_value};
      }
      Iterator& operator++() /* prefix increment*/
      {
        next_char_ptr++;
        last_qgram_was_processed = false;
        return *this;
      }
      bool operator != (const Iterator& other) const
      {
        return static_cast<bool>(next_char_ptr < other.end_of_sequence);
      }
  };

    QgramTransformer qgram_transformer;
    size_t qgram_length;
    const SequenceBaseType *sequence;
    size_t seqlen;
    CyclicBuffer_uint8 current_window;
    uint64_t max_integer_code;
#ifndef NDEBUG
    uint8_t qgram_buffer[MAX_QGRAM_LENGTH];
#endif
  public:
    QgramRecHashValueIterator(size_t _qgram_length,
                              const SequenceBaseType *_sequence,
                              size_t _seqlen)
      : qgram_transformer(QgramTransformer(_qgram_length))
      , qgram_length(_qgram_length)
      , sequence(_sequence)
      , seqlen(_seqlen)
    {
      assert(qgram_length <= 32);
      max_integer_code = qgram_length == 32
                           ? UINT64_MAX
                           : std::pow(alpha_size,_qgram_length) - 1;
      current_window.initialize(_qgram_length);
    }
    Iterator begin(void)
    {
      uint64_t this_hash_value;
      uint64_t this_compl_hash_value;
      if (qgram_length <= seqlen)
      {
        for (const SequenceBaseType *qgram_ptr = sequence + qgram_length - 1;
             qgram_ptr >= sequence; qgram_ptr--)
        {
          const uint8_t rank = alphabet.char_to_rank(*qgram_ptr);
          current_window.prepend(rank);
        }
        std::tie(this_hash_value,
                 this_compl_hash_value)
          = qgram_transformer.first_hash_value_pair_get(
                                 current_window.pointer_to_array(),
                                 qgram_length);
      } else
      {
        /* the next two values will not be used as there is no qgram */
        this_hash_value = this_compl_hash_value = 0;
      }
      return Iterator(qgram_transformer,
                      qgram_length,
                      current_window,
                      sequence,
                      this_hash_value,
                      this_compl_hash_value);
    }
    Iterator end(void)
    {
      return Iterator(qgram_transformer,
                      current_window,
                      sequence + seqlen);
    }
    /* The following functions are for qgram invertible integer codes only */
    uint64_t qgram_encode(const SequenceBaseType *qgram)
    {
      uint64_t code = 0;
      uint64_t mult= 1;
      for (const SequenceBaseType *qgram_ptr = qgram + qgram_length - 1;
           qgram_ptr >= qgram; qgram_ptr--)
      {
        const uint8_t rank = alphabet.char_to_rank(*qgram_ptr);
        current_window.prepend(rank);
        code += mult * static_cast<uint64_t>(rank);
        mult *= alpha_size;
      }
      return code;
    }
    uint64_t qgram_encode_complement(const SequenceBaseType *qgram)
    {
      uint64_t code = 0;
      uint64_t mult= 1;
      for (const SequenceBaseType *qgram_ptr = qgram;
           qgram_ptr < qgram + qgram_length; qgram_ptr++)
      {
        const uint8_t rank = alphabet.char_to_rank(*qgram_ptr);
        code += mult * static_cast<uint64_t>(rank);
        mult *= alpha_size;
      }
      return code;
    }
    [[nodiscard]] size_t qgram_length_get(void) const noexcept
    {
      return qgram_length;
    }
#ifndef NDEBUG
    const uint8_t *qgram_decode(uint64_t code)
    {
      assert(code <= max_integer_code);
      for (uint8_t *qgram_ptr = qgram_buffer + qgram_length - 1;
           qgram_ptr >= qgram_buffer;
           qgram_ptr--)
      {
        *qgram_ptr = code % alpha_size;
        code /= alpha_size;
      }
      return &qgram_buffer[0];
    }
    void qgrams_compare(const SequenceBaseType *orig_qgram,
                        const uint8_t *mapped_qgram)
    {
      for (size_t idx = 0; idx < qgram_length; idx++)
      {
        uint8_t rank = alphabet.char_to_rank(orig_qgram[idx]);
        assert(rank == mapped_qgram[idx]);
      }
    }
#endif
};
#endif
