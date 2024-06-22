/*
  Copyright (c) 2021 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2021 Center for Bioinformatics, University of Hamburg

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

#ifndef DNA_SEQ_ENCODER_HPP
#define DNA_SEQ_ENCODER_HPP

#include <cassert>
#include <cstdlib>
#include <cstdint>
#include <cinttypes>
#include <algorithm>
#ifndef NDEBUG
#include <iostream>
#endif
#include "utilities/mathsupport.hpp"
#include "sequences/gttl_fastq_iterator.hpp"
#include "sequences/alphabet.hpp"

/* Wildcards are transformed to rank 0 */
static constexpr const alphabet::GttlAlphabet_UL_0 dna_alphabet;

template<typename StoreUnitType>
class DNASeqEncoder
{
  private:
  /* Number of bits which can be stored in a unit */
  static constexpr const int bits_in_store_unit
    = sizeof(StoreUnitType) * CHAR_BIT;
  /* Number of character which can be stored in a unit:
     we need two bits per character */
  static constexpr const int characters_per_unit = bits_in_store_unit/2;
  const size_t prefix_length;
  const int additional_bits;
#ifndef NDEBUG
  const size_t additional_value_max;
#endif
  const size_t num_bits;
  const size_t num_units;
  const int additional_shift;
  const StoreUnitType end_mask;
  StoreUnitType end_mask_get(void) const noexcept
  {
    const size_t endbits = bits_in_store_unit -
                           2 * (prefix_length % characters_per_unit);
    assert(endbits > 0);
    return (static_cast<StoreUnitType>(1) << endbits) - 1;
  }
  int additional_shift_get(void) const noexcept
  {
    const int remainder = prefix_length % characters_per_unit;
    if (remainder > 0)
    {
      if (additional_bits >= (bits_in_store_unit - remainder * 2))
      {
        return additional_bits - (bits_in_store_unit - remainder * 2);
      } else
      {
        return 0;
      }
    } else
    {
      return additional_bits;
    }
  }
  public:
  DNASeqEncoder(size_t _prefix_length,int _additional_bits = 0)
    : prefix_length(_prefix_length)
    , additional_bits(_additional_bits)
#ifndef NDEBUG
    , additional_value_max((size_t(1) << _additional_bits) - 1)
#endif
    , num_bits(std::max(size_t(64),2 * prefix_length + additional_bits))
    , num_units((num_bits + bits_in_store_unit - 1)/bits_in_store_unit)
    , additional_shift(additional_shift_get())
    , end_mask(end_mask_get())
  {
#ifndef NDEBUG
    std::cout << "# num_bits=" << num_bits << std::endl;
    std::cout << "# num_units=" << num_units << std::endl;
    assert(additional_bits
             < static_cast<int>(sizeof(size_t) * bits_in_store_unit));
    std::cout << "# additional_bits=" << additional_bits << std::endl;
    std::cout << "# bits_in_store_unit=" << bits_in_store_unit << std::endl;
    std::cout << "# prefix_length=" << prefix_length << std::endl;
    std::cout << "# characters_per_unit=" << characters_per_unit << std::endl;
#endif
  }
  void encode(StoreUnitType *encoding, const char *char_seq, size_t a_val = 0)
       const noexcept
  {
    assert(a_val <= additional_value_max);
    size_t encoding_index = 0, char_idx = 0;

    if (prefix_length >= characters_per_unit)
    {
      while (char_idx + characters_per_unit <= prefix_length)
      {
        int shift = bits_in_store_unit - 2;
        StoreUnitType value = 0;
        for (size_t idx = 0; idx < characters_per_unit; idx++)
        {
          const uint8_t r = dna_alphabet.char_to_rank(char_seq[char_idx+idx]);
          value |= (static_cast<StoreUnitType>(r)) << shift;
          shift -= 2;
        }
        assert(encoding_index < num_units);
        encoding[encoding_index++] = value;
        char_idx += characters_per_unit;
      }
    }
    assert(prefix_length < char_idx + characters_per_unit);
    if (prefix_length % characters_per_unit > 0)
    {
      int shift = bits_in_store_unit - 2;
      StoreUnitType value = 0;
      for(size_t idx = 0; idx < prefix_length % characters_per_unit; idx++)
      {
        const uint8_t r = dna_alphabet.char_to_rank(char_seq[char_idx+idx]);
        value |= (static_cast<StoreUnitType>(r)) << shift;
        shift -= 2;
      }
      assert(encoding_index < num_units);
      encoding[encoding_index++] = value;
    }
    int remaining_additional_bits = additional_shift;
    while (remaining_additional_bits >= bits_in_store_unit)
    {
      remaining_additional_bits -= bits_in_store_unit;
      assert(encoding_index < num_units);
      encoding[encoding_index++]
        = static_cast<StoreUnitType>(a_val >> remaining_additional_bits);
    }
    if (remaining_additional_bits > 0)
    {
      assert(remaining_additional_bits < bits_in_store_unit);
      encoding[encoding_index++]
        = static_cast<StoreUnitType>(a_val << (bits_in_store_unit -
                                               remaining_additional_bits));
    }
  }
#ifndef NDEBUG
  void sequence_encoding_verify(const StoreUnitType *encoding,
                                const char *char_seq)
       const noexcept
  {
    size_t encoding_index = 0;
    int shift = bits_in_store_unit - 2;
    for (size_t char_idx = 0; char_idx < prefix_length; char_idx++)
    {
      const uint8_t r = dna_alphabet.char_to_rank(char_seq[char_idx]),
                   er = static_cast<uint8_t>((encoding[encoding_index] >> shift)
                                             & uint8_t(3));
      if (r != er)
      {
        std::cerr << "r = " << static_cast<int>(r) << " != "
                  << static_cast<int>(er) << " = er" << std::endl;
        exit(EXIT_FAILURE);
      }
      if (shift >= 2)
      {
        shift -= 2;
      } else
      {
        assert(shift == 0);
        shift = bits_in_store_unit - 2;
        encoding_index++;
      }
    }
  }
  std::string decode(const StoreUnitType *encoding) const noexcept
  {
    std::string s;
    for (size_t idx = 0; idx < num_units; idx++)
    {
      for (int shift = bits_in_store_unit - 8; shift >= 0; shift -= 8)
      {
        s += std::to_string((encoding[idx] >> shift) &
                            static_cast<StoreUnitType>(UINT8_MAX));
        if (shift > 0)
        {
          s += std::string(" ");
        }
      }
      if (idx < num_units-1)
      {
        s += std::string(" ");
      }
    }
    return s;
  }
#endif
  size_t decode_additional_value(const StoreUnitType *encoding) const noexcept
  {
    size_t additional_value;
    if (additional_shift < additional_bits)
    {
      additional_value = static_cast<size_t>(encoding[prefix_length/4] &
                                             end_mask)
                         << additional_shift;
    } else
    {
      additional_value = 0;
    }
    int remaining_bits = additional_shift;
    size_t decoding_index = (prefix_length+3)/4;
    while (remaining_bits >= bits_in_store_unit)
    {
      remaining_bits -= bits_in_store_unit;
      assert(decoding_index < num_units);
      additional_value |= static_cast<size_t>(encoding[decoding_index++]
                                              << remaining_bits);
    }
    if (remaining_bits > 0)
    {
      assert(remaining_bits < bits_in_store_unit && decoding_index < num_units);
      additional_value |= static_cast<size_t>(encoding[decoding_index]
                                              >> (bits_in_store_unit -
                                                  remaining_bits));
    }
    return additional_value;
  }
  size_t num_bits_get(void) const noexcept
  {
    return num_bits;
  }
  size_t num_units_get(void) const noexcept
  {
    return num_units;
  }
  size_t num_sequence_bytes_get(void) const noexcept
  {
    const size_t num_sequence_bits(std::max(size_t(64),2 * prefix_length));
    return (num_sequence_bits + bits_in_store_unit - 1)/bits_in_store_unit;
  }
};

template<typename StoreUnitType>
class DNAEncodingForLength
{
  const DNASeqEncoder<StoreUnitType> dna_seq_encoder;
  const size_t constant_sequence_length, num_units;
  size_t allocated, nextfree, add_factor;
  StoreUnitType *units;
  StoreUnitType *append_ptr(void)
  {
    if (nextfree + num_units >= allocated)
    {
      allocated += num_units * add_factor;
      units = static_cast<StoreUnitType *>
                         (realloc(units,allocated * sizeof *units));
      add_factor *= 1.8;
    }
    StoreUnitType *ptr = units + nextfree;
    nextfree += num_units;
    return ptr;
  }
  public:
  DNAEncodingForLength(size_t _constant_sequence_length)
    : dna_seq_encoder(_constant_sequence_length)
    , constant_sequence_length(_constant_sequence_length)
    , num_units(dna_seq_encoder.num_units_get())
    , allocated(0)
    , nextfree(0)
    , add_factor(1000)
    , units(nullptr)
  {
  }
  ~DNAEncodingForLength(void)
  {
    free(units);
  }
  void add(const std::string_view &sequence)
  {
    StoreUnitType *ptr = append_ptr();
    dna_seq_encoder.encode(ptr,sequence.data());
#ifndef NDEBUG
    dna_seq_encoder.sequence_encoding_verify(ptr,sequence.data());
#endif
  }
  void final_resize(void)
  {
    units = static_cast<StoreUnitType *>
                       (realloc(units,nextfree * sizeof *units));
    allocated = nextfree;
  }
  size_t num_units_get(void) const
  {
    return num_units;
  }
  size_t number_of_sequences_get(void) const
  {
    assert(nextfree % num_units == 0);
    return nextfree / num_units;
  }
  size_t sequence_length_get(void) const
  {
    return constant_sequence_length;
  }
  const StoreUnitType *units_get(void) const
  {
    return units;
  }
  void statistics(void) const
  {
    std::cout << "# length of sequences\t"
              << sequence_length_get() << std::endl;
    std::cout << "# number of sequences\t"
              << number_of_sequences_get() << std::endl;
    std::cout << "# units per sequence\t"
              << num_units_get() << std::endl;
    std::cout << "# total size (MB)\t"
              << static_cast<size_t>(mega_bytes(number_of_sequences_get() *
                                                num_units_get() *
                                                sizeof(StoreUnitType)))
              << std::endl;
  }
  std::string to_string(void) const
  {
    static const std::array<char,4> dna_letters{'A','C','G','T'};
    static constexpr const int bits_in_store_unit
      = sizeof(StoreUnitType) * CHAR_BIT;
    int shift = bits_in_store_unit - 2;
    std::string s;
    size_t unit_num = 0;
    for (size_t idx = 0; idx < this->sequence_length_get(); idx++)
    {
      const size_t char_idx = static_cast<size_t>(units[unit_num] >> shift)
                              & static_cast<size_t>(3);
      s += dna_letters[char_idx];
      if (shift > 0)
      {
        assert(shift > 1);
        shift -= 2;
      } else
      {
        unit_num++;
        shift = bits_in_store_unit - 2;
      }
    }
    return s;
  }
};

template<typename StoreUnitType>
class DNAEncodingMultiLength
{
  using ThisDNAEncodingForLength = DNAEncodingForLength<StoreUnitType>;
  std::vector<ThisDNAEncodingForLength *> enc_vec;
  public:
  DNAEncodingMultiLength(const std::string &inputfilename)
  {
    constexpr const int buf_size = 1 << 14;
    GttlLineIterator<buf_size> line_iterator(inputfilename.c_str());
    GttlFastQIterator<GttlLineIterator<buf_size>> fastq_it(line_iterator);
    for (auto &fastq_entry : fastq_it)
    {
      const std::string_view &sequence = fastq_entry.sequence_get();
      for (size_t idx = enc_vec.size(); idx <= sequence.size(); idx++)
      {
        enc_vec.push_back(nullptr);
      }
      ThisDNAEncodingForLength *enc_ptr;
      assert(sequence.size() < enc_vec.size());
      if (enc_vec[sequence.size()] != nullptr)
      {
        enc_ptr = enc_vec[sequence.size()];
      } else
      {
        enc_ptr = new ThisDNAEncodingForLength(sequence.size());
      }
      enc_ptr->add(sequence);
      enc_vec[sequence.size()] = enc_ptr;
    }
    size_t w_idx = 0;
    for (size_t idx = 0; idx < enc_vec.size(); idx++)
    {
      if (enc_vec[idx] != nullptr)
      {
        enc_vec[idx]->final_resize();
        assert(w_idx < idx);
        enc_vec[w_idx++] = enc_vec[idx];
      }
    }
    assert(w_idx > 0);
    enc_vec.resize(w_idx);
  }
  ~DNAEncodingMultiLength(void)
  {
    for (auto &dna_encoding : enc_vec)
    {
      assert(dna_encoding != nullptr);
      delete dna_encoding;
    }
  }
  void statistics(void) const
  {
    for (auto &dna_encoding : enc_vec)
    {
      assert(dna_encoding != nullptr);
      dna_encoding->statistics();
    }
  }
  class Iterator
  {
    const std::vector<ThisDNAEncodingForLength *> &enc_vec_ref;
    size_t current_enc_vec_idx,
           current_seqnum;
    bool exhausted;
    public:
    Iterator(const std::vector<ThisDNAEncodingForLength *> &_enc_vec_ref,
             bool _exhausted)
      : enc_vec_ref(_enc_vec_ref)
      , current_enc_vec_idx(0)
      , current_seqnum(0)
      , exhausted(_exhausted)
    { }
    std::pair<const uint64_t *,size_t> operator *(void) const
    {
      const uint64_t *units
        = enc_vec_ref[current_enc_vec_idx]->units_get();
      const size_t num_units = enc_vec_ref[current_enc_vec_idx]
                                 ->num_units_get();
      return std::make_pair(units + current_seqnum * num_units,
                            enc_vec_ref[current_enc_vec_idx]
                              ->sequence_length_get());
    }
    bool operator != (const Iterator& other) const noexcept
    {
      return exhausted != other.exhausted;
    }
    Iterator& operator++() /* prefix increment*/
    {
      assert(not exhausted);
      if (current_seqnum + 1 <
          enc_vec_ref[current_enc_vec_idx]->number_of_sequences_get())
      {
        current_seqnum++;
      } else
      {
        current_seqnum = 0;
        if (current_enc_vec_idx + 1 < enc_vec_ref.size())
        {
          current_enc_vec_idx++;
        } else
        {
          exhausted = true;
        }
      }
      return *this;
    }
  };
  Iterator begin() const
  {
    return Iterator(enc_vec,false);
  }
  Iterator end() const
  {
    return Iterator(enc_vec,true);
  }
};

template<typename StoreUnitType>
class DNAEncoding
{
  size_t constant_sequence_length, num_units, allocated, nextfree,
         add_factor;
  StoreUnitType *units;
  StoreUnitType *append_ptr(void)
  {
    if (nextfree + num_units >= allocated)
    {
      allocated += num_units * add_factor;
      units = static_cast<StoreUnitType *>
                         (realloc(units,allocated * sizeof *units));
      add_factor *= 1.1;
    }
    StoreUnitType *ptr = units + nextfree;
    nextfree += num_units;
    return ptr;
  }
  public:
  DNAEncoding(const std::string &inputfilename)
    : constant_sequence_length(0)
    , num_units(0)
    , allocated(0)
    , nextfree(0)
    , add_factor(1000000)
    , units(nullptr)
  {
    constexpr const int buf_size = 1 << 14;
    GttlLineIterator<buf_size> line_iterator(inputfilename.c_str());
    GttlFastQIterator<GttlLineIterator<buf_size>> fastq_it(line_iterator);
    auto fastq_entry = fastq_it.begin();
    const std::string_view &first_sequence = (*fastq_entry).sequence_get();
    constant_sequence_length = first_sequence.size();
    DNASeqEncoder<StoreUnitType> dna_seq_encoder(constant_sequence_length);
    num_units = dna_seq_encoder.num_units_get();
    StoreUnitType *ptr = append_ptr();
    dna_seq_encoder.encode(ptr,first_sequence.data());
#ifndef NDEBUG
    dna_seq_encoder.sequence_encoding_verify(ptr,first_sequence.data());
#endif
    ++fastq_entry;
    while (fastq_entry != fastq_it.end())
    {
      const std::string_view &sequence = (*fastq_entry).sequence_get();
      if (sequence.size() != constant_sequence_length)
      {
        throw std::string("can only handle sequence of same length: "
                          "first sequence has length ") +
              std::to_string(constant_sequence_length) +
              std::string(" while this sequence has length ") +
              std::to_string(sequence.size());
      }
      ptr = append_ptr();
      dna_seq_encoder.encode(ptr,sequence.data());
#ifndef NDEBUG
      dna_seq_encoder.sequence_encoding_verify(ptr,sequence.data());
#endif
      ++fastq_entry;
    }
    units = static_cast<StoreUnitType *>
                       (realloc(units,nextfree * sizeof *units));
    allocated = nextfree;
  }
  ~DNAEncoding(void)
  {
    free(units);
  }
  size_t num_units_get(void) const
  {
    return num_units;
  }
  size_t number_of_sequences_get(void) const
  {
    assert(nextfree % num_units == 0);
    return nextfree / num_units;
  }
  size_t sequence_length_get(void) const
  {
    return constant_sequence_length;
  }
  const StoreUnitType *units_get(void) const
  {
    return units;
  }
  void statistics(void) const
  {
    std::cout << "# length of sequences\t"
              << sequence_length_get() << std::endl;
    std::cout << "# number of sequences\t"
              << number_of_sequences_get() << std::endl;
    std::cout << "# units per sequence\t"
              << num_units_get() << std::endl;
    std::cout << "# total size (MB)\t"
              << static_cast<size_t>(mega_bytes(number_of_sequences_get() *
                                                num_units_get() *
                                                sizeof(StoreUnitType)))
              << std::endl;
  }
  std::string to_string(void) const
  {
    static const std::array<char,4> dna_letters{'A','C','G','T'};
    static constexpr const int bits_in_store_unit
      = sizeof(StoreUnitType) * CHAR_BIT;
    int shift = bits_in_store_unit - 2;
    std::string s;
    size_t unit_num = 0;
    for (size_t idx = 0; idx < this->sequence_length_get(); idx++)
    {
      const size_t char_idx = static_cast<size_t>(units[unit_num] >> shift)
                              & static_cast<size_t>(3);
      s += dna_letters[char_idx];
      if (shift > 0)
      {
        assert(shift > 1);
        shift -= 2;
      } else
      {
        unit_num++;
        shift = bits_in_store_unit - 2;
      }
    }
    return s;
  }
};

class DNAQgramDecoder
{
  class Iterator
  {
    const uint64_t mask;
    const uint64_t *sub_unit_ptr;
    size_t current_qgram_idx, idx_of_unit;
    size_t shift_last;
    uint64_t integer;
    public:
    Iterator(size_t qgram_length,
             size_t _current_qgram_idx,
             const uint64_t *_sub_unit_ptr)
      : mask(((uint64_t(1) << (2 * (qgram_length - 1))) - 1) << 2)
      , sub_unit_ptr(_sub_unit_ptr)
      , current_qgram_idx(_current_qgram_idx)
      , idx_of_unit(_current_qgram_idx + qgram_length < size_t(32) ? 0
                                                                   : size_t(1))
      , shift_last(qgram_length == size_t(32) ? size_t(62)
                                              : (size_t(62) -
                                                 (2 * qgram_length)))
      , integer(sub_unit_ptr[0] >> (size_t(64) - 2 * qgram_length))
    { }
    uint64_t operator*() const
    {
      return integer;
    }
    Iterator& operator++() /* prefix increment*/
    {
      integer = (integer << 2) & mask;
      assert((integer & uint64_t(3)) == 0);
      const uint64_t new_char = (sub_unit_ptr[idx_of_unit] >> shift_last)
                                & uint64_t(3);
      integer |= new_char;
      if (shift_last > 0)
      {
        assert(shift_last >= 2);
        shift_last -= 2;
      } else
      {
        shift_last = 62;
        idx_of_unit++;
      }
      current_qgram_idx++;
      return *this;
    }
    bool operator != (const Iterator& other) const
    {
      return current_qgram_idx != other.current_qgram_idx;
    }
  };
  const size_t qgram_length;
  const size_t number_of_qgrams;
  const uint64_t *sub_unit_ptr;
  public:
  DNAQgramDecoder(size_t _qgram_length,
                  size_t _number_of_qgrams,
                  const uint64_t *_sub_unit_ptr)
    : qgram_length(_qgram_length)
    , number_of_qgrams(_number_of_qgrams)
    , sub_unit_ptr(_sub_unit_ptr)
  {}
  Iterator begin(void) const
  {
    return Iterator(qgram_length,0,sub_unit_ptr);
  }
  Iterator end(void) const
  {
    return Iterator(qgram_length,number_of_qgrams,sub_unit_ptr);
  }
};
#endif
