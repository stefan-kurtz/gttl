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
#ifndef NDEBUG
#include <iostream>
#endif
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
    std::cout << "num_bits=" << num_bits << std::endl;
    std::cout << "num_units=" << num_units << std::endl;
    assert(additional_bits
             < static_cast<int>(sizeof(size_t) * bits_in_store_unit));
    std::cout << "additional_bits=" << additional_bits << std::endl;
    std::cout << "bits_in_store_unit=" << bits_in_store_unit << std::endl;
    std::cout << "prefix_length=" << prefix_length << std::endl;
    std::cout << "characters_per_unit=" << characters_per_unit << std::endl;
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
};

template<typename StoreUnitType>
class DNAEncoding
{
  size_t constant_sequence_length, num_units, allocated, nextfree;
  StoreUnitType *units;
  StoreUnitType *append_ptr(size_t factor)
  {
    if (nextfree + num_units >= allocated)
    {
      allocated += num_units * factor;
      units = static_cast<StoreUnitType *>
                         (realloc(units,allocated * sizeof *units));
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
    const size_t space_factor = 1000;
    StoreUnitType *ptr = append_ptr(space_factor);
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
      ptr = append_ptr(space_factor);
      dna_seq_encoder.encode(ptr,sequence.data());
#ifndef NDEBUG
      dna_seq_encoder.sequence_encoding_verify(ptr,sequence.data());
#endif
      ++fastq_entry;
    }
  }
  ~DNAEncoding(void)
  {
    free(units);
  }
  void reset(void)
  {
    free(units);
    units = nullptr;
    nextfree = allocated = num_units = 0;
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
  size_t sizeof_unit_get() const
  {
    return sizeof(StoreUnitType);
  }
  void statistics(void) const
  {
    std::cout << "# number of sequences\t"
              << number_of_sequences_get() << std::endl;
    std::cout << "# length of sequences\t"
              << sequence_length_get() << std::endl;
    std::cout << "# units per sequence\t"
              << num_units_get() << std::endl;
    std::cout << "# total size (MB)\t"
              << static_cast<size_t>(mega_bytes(number_of_sequences_get() *
                                                num_units_get() *
                                                sizeof_unit_get()))
              << std::endl;
  }
};
#endif
