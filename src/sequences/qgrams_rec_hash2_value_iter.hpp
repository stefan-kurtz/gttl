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
#ifndef QGRAMS_HASH2_REC_ITER_HPP
#define QGRAMS_HASH2_REC_ITER_HPP
#include <cstdbool>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include "utilities/cyclic_buffer.hpp"
#include "sequences/max_qgram_length.hpp"

static uint8_t nucleotide2rank(char cc)
{
  return (static_cast<uint8_t>(cc) >> 1) & uint8_t(3);
}

/* Implementation of iterator class follows concept described in
   https://davidgorski.ca/posts/stl-iterators/ */

template<uint64_t first_hash_value_get(const uint8_t *,size_t),
         uint64_t first_compl_hash_value_get(const uint8_t *,size_t),
         class AuxData,
         AuxData aux_data_get(size_t),
         uint64_t next_hash_value_get(uint8_t,uint64_t,uint8_t,AuxData &),
         uint64_t next_compl_hash_value_get(uint8_t,uint64_t,uint8_t,AuxData &)>
class QgramRecHash2ValueIterator
{
  static constexpr const size_t alpha_size = 4;
  using SequenceBaseType = char;
  using CyclicBuffer_uint8 = CyclicBuffer<uint8_t,MAX_QGRAM_LENGTH>;

  struct Iterator
  {
    private:
      CyclicBuffer_uint8 &current_window;
      const SequenceBaseType *next_char_ptr;
      bool last_qgram_was_processed;
      const SequenceBaseType *end_of_sequence;
      AuxData aux_data;
      uint64_t hash_value, compl_hash_value;
      uint8_t wildcards_in_qgram;
    public:
      using iterator_category = std::forward_iterator_tag;
      using value_type = SequenceBaseType;
      using difference_type = size_t;
      using pointer = const SequenceBaseType*;
      using reference = uint64_t&;

      /* Constructor for begin() */
      Iterator(size_t _qgram_length,
               CyclicBuffer_uint8 &_current_window,
               pointer _sequence,
               uint64_t first_hash_value,
               uint64_t first_compl_hash_value)
        : current_window(_current_window)
        , next_char_ptr(_sequence + _qgram_length - 1)
        , last_qgram_was_processed(true)
        , end_of_sequence(nullptr)
        , aux_data(aux_data_get(_qgram_length))
        , hash_value(first_hash_value)
        , compl_hash_value(first_compl_hash_value)
      {
      }
      /* Constructor for end() */
      Iterator(CyclicBuffer_uint8 &_current_window,
               pointer _end_of_sequence) :
        current_window(_current_window),
        end_of_sequence(_end_of_sequence)
      {
      }
      std::pair<uint64_t,uint64_t> operator*()
      {
        if (!last_qgram_was_processed)
        {
          const uint8_t new_rank = nucleotide2rank(*next_char_ptr);
          const uint8_t old_rank = current_window.shift(new_rank);
          hash_value = next_hash_value_get(old_rank,
                                           hash_value,
                                           new_rank,
                                           aux_data);
          compl_hash_value = next_compl_hash_value_get(old_rank,
                                                       hash_value,
                                                       new_rank,
                                                       aux_data);
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

  private:
    size_t qgram_length;
    const SequenceBaseType *sequence;
    size_t seqlen;
    CyclicBuffer_uint8 current_window;
    uint64_t max_integer_code;
#ifndef NDEBUG
    uint8_t qgram_buffer[MAX_QGRAM_LENGTH];
#endif
  public:
    QgramRecHash2ValueIterator(size_t _qgram_length,
                               const SequenceBaseType *_sequence,
                               size_t _seqlen):
      qgram_length(_qgram_length),
      sequence(_sequence),
      seqlen(_seqlen)
    {
      max_integer_code = qgram_length == 32
                           ? UINT64_MAX
                           : std::pow(alpha_size,_qgram_length) - 1;
      current_window.initialize(_qgram_length);
    }
    Iterator begin()
    {
      uint64_t this_hash_value, this_compl_hash_value;
      if (qgram_length <= seqlen)
      {
        for (const SequenceBaseType *qgram_ptr = sequence + qgram_length - 1;
             qgram_ptr >= sequence; qgram_ptr--)
        {
          const uint8_t rank = nucleotide2rank(*qgram_ptr);
          current_window.prepend(rank);
        }
        this_hash_value
          = first_hash_value_get(current_window.pointer_to_array(),
                                 qgram_length);
        this_compl_hash_value
          = first_compl_hash_value_get(current_window.pointer_to_array(),
                                       qgram_length);
      } else
      {
        /* the next two values will not be used as there is no qgram */
        this_hash_value = this_compl_hash_value = 0;
      }
      return Iterator(qgram_length,current_window,sequence,
                      this_hash_value,this_compl_hash_value);
    }
    Iterator end()
    {
      return Iterator(current_window,sequence + seqlen);
    }
    /* The following functions are for qgram integer codes only */
    uint64_t qgram_encode(const SequenceBaseType *qgram)
    {
      uint64_t code = 0, mult= 1;
      for (const SequenceBaseType *qgram_ptr = qgram + qgram_length - 1;
           qgram_ptr >= qgram; qgram_ptr--)
      {
        const uint8_t rank = nucleotide2rank(*qgram_ptr);
        current_window.prepend(rank);
        code += mult * static_cast<uint64_t>(rank);
        mult *= alpha_size;
      }
      return code;
    }
    uint64_t qgram_encode_complement(const SequenceBaseType *qgram)
    {
      uint64_t code = 0, mult= 1;
      for (const SequenceBaseType *qgram_ptr = qgram;
           qgram_ptr < qgram + qgram_length; qgram_ptr++)
      {
        const uint8_t rank = nucleotide2rank(*qgram_ptr);
        code += mult * static_cast<uint64_t>(rank);
        mult *= alpha_size;
      }
      return code;
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
        uint8_t rank = nucleotide2rank(orig_qgram[idx]);
        assert(rank == mapped_qgram[idx]);
      }
    }
#endif
};
#endif