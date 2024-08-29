/*
  Copyright (c) 2024 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2024 Center for Bioinformatics, University of Hamburg

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
#include <map>
#ifndef NDEBUG
#include <iostream>
#endif
#include "utilities/mathsupport.hpp"
#include "utilities/unused.hpp"
#include "sequences/gttl_fastq_iterator.hpp"
#include "sequences/alphabet.hpp"

template<typename StoreUnitType,bool verbose>
class DNASeqEncoder
{
  private:
  /* Wildcards are transformed to rank 0 */
  static constexpr const alphabet::GttlAlphabet_UL_0 dna_alphabet{};
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
    if constexpr (verbose)
    {
      std::cout << "# num_bits=" << num_bits << std::endl;
      std::cout << "# num_units=" << num_units << std::endl;
      assert(additional_bits
               < static_cast<int>(sizeof(size_t) * bits_in_store_unit));
      std::cout << "# additional_bits=" << additional_bits << std::endl;
      std::cout << "# bits_in_store_unit=" << bits_in_store_unit << std::endl;
      std::cout << "# prefix_length=" << prefix_length << std::endl;
      std::cout << "# characters_per_unit=" << characters_per_unit << std::endl;
    }
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
    const size_t decoding_index = (prefix_length+3)/4;
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

template<typename StoreUnitType,bool verbose>
class DNAEncodingForLength
{
  const DNASeqEncoder<StoreUnitType,verbose> dna_seq_encoder;
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
  size_t num_units_get(void) const noexcept
  {
    return num_units;
  }
  size_t number_of_sequences_get(void) const noexcept
  {
    assert(nextfree % num_units == 0);
    return nextfree / num_units;
  }
  size_t sequence_length_get(void) const noexcept
  {
    return constant_sequence_length;
  }
  const StoreUnitType *units_get(void) const noexcept
  {
    return units;
  }
  size_t total_size_get(void) const noexcept
  {
    return number_of_sequences_get() * num_units_get();
  }
  void statistics(void) const noexcept
  {
    std::cout << "# length of sequences\t"
              << sequence_length_get() << std::endl;
    std::cout << "# number of sequences\t"
              << number_of_sequences_get() << std::endl;
    std::cout << "# units per sequence\t"
              << num_units_get() << std::endl;
    std::cout << "# total size (MB)\t"
              << static_cast<size_t>(mega_bytes(total_size_get() *
                                                sizeof(StoreUnitType)))
              << std::endl;
  }
  std::string to_string(void) const noexcept
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

static bool decide_append_previous(const std::vector<size_t> &size_vec,
                                   size_t total_size,
                                   size_t local_sum)
{
  assert(size_vec.size() > 0);
  const double mean0 = static_cast<double>(total_size)/size_vec.size();
  double sum_squared_difference = 0;
  for (size_t idx = 0; idx < size_vec.size(); idx++)
  {
    const size_t s = size_vec[idx];
    const double diff = (idx < size_vec.size() - 1 ? s
                                                   : (s + local_sum)) - mean0;
    sum_squared_difference += (diff * diff);
  }
  const double var0 = sum_squared_difference/size_vec.size();
  const double mean1 = total_size/(size_vec.size() + 1);
  sum_squared_difference = 0;
  for (auto s : size_vec)
  {
    const double diff = s - mean1;
    sum_squared_difference += (diff * diff);
  }
  const double diff = local_sum - mean1;
  sum_squared_difference += diff * diff;
  const double var1 = sum_squared_difference/(size_vec.size() + 1);
  return var0 < var1;
}


template<typename StoreUnitType,bool verbose>
class DNAEncodingMultiLength
{
  using ThisDNAEncodingForLength = DNAEncodingForLength<StoreUnitType,verbose>;
  std::vector<ThisDNAEncodingForLength *> enc_vec;
  size_t total_size, total_number_of_sequences;
  using KeyValuesType = std::vector<std::tuple<size_t,size_t,size_t,size_t>>;
  KeyValuesType expanded_vec;
  std::vector<size_t> end_idx_of_part_vec;
  auto key_values_vector_get(void) const
  {
    KeyValuesType key_values;
    size_t enc_vec_idx = 0;
    for (auto &dna_encoding : enc_vec)
    {
      key_values.push_back(std::make_tuple(enc_vec_idx,
                                           dna_encoding
                                             ->number_of_sequences_get(),
                                           dna_encoding->num_units_get(),
                                           0));
      enc_vec_idx++;
    }
    return key_values;
  }
  void key_values_show(const KeyValuesType &key_values) const
  {
    size_t idx = 0;
    for (auto &kv : key_values)
    {
      std::cout << idx << "\t"
                << enc_vec[std::get<0>(kv)]->sequence_length_get() << "\t"
                << std::get<1>(kv) << "\t"
                << std::get<2>(kv) << "\t"
                << std::get<3>(kv) << "\t"
                << (std::get<1>(kv) * std::get<2>(kv))
                << std::endl;
      idx++;
    }
  }
  void divide_vector_evenly(size_t num_parts)
  {
    assert(num_parts > 1);
    const size_t mean = (total_size + num_parts - 1)/num_parts;
    if constexpr (verbose)
    {
      std::cout << "# size_sum\t" << total_size << std::endl;
      std::cout << "# num_parts\t" << num_parts << std::endl;
      std::cout << "# mean " << mean << std::endl;
    }
    size_t local_sum = 0, idx = 0;
    std::vector<size_t> size_vec;
    assert(end_idx_of_part_vec.size() == 0);
    for (auto &kv : expanded_vec)
    {
      auto s = std::get<1>(kv) * std::get<2>(kv);
      if (local_sum + s <= mean)
      {
        local_sum += s;
      } else
      {
        if (local_sum > 0)
        {
          end_idx_of_part_vec.push_back(idx);
          size_vec.push_back(local_sum);
        }
        local_sum = s;
      }
      idx++;
    }
    if (local_sum > 0)
    {
      if (decide_append_previous(size_vec,total_size,local_sum))
      {
        end_idx_of_part_vec.back() = expanded_vec.size();
        size_vec.back() += local_sum;
      } else
      {
        end_idx_of_part_vec.push_back(expanded_vec.size());
        size_vec.push_back(local_sum);
      }
    }
    size_t local_sum_sum = 0, begin_idx = 0;
    assert(size_vec.size() == end_idx_of_part_vec.size());
    for (size_t j = 0; j < size_vec.size(); j++)
    {
      const size_t end_idx = end_idx_of_part_vec.at(j);
      const size_t local_sum = size_vec[j];
      if constexpr (verbose)
      {
        std::cout << end_idx << "\t" << local_sum << std::endl;
      }
      size_t this_sum = 0;
      for (size_t idx = begin_idx; idx < end_idx; idx++)
      {
        this_sum += std::get<1>(expanded_vec[idx]) *
                    std::get<2>(expanded_vec[idx]);
      }
      if (local_sum != this_sum)
      {
        std::cerr << "local_sum = " << local_sum << " != " << this_sum
                  << " = this_sum" << std::endl;
        exit(EXIT_FAILURE);
      }
      local_sum_sum += local_sum;
      begin_idx = end_idx;
    }
    if (local_sum_sum != total_size)
    {
      std::cerr << "local_sum_sum = " << local_sum_sum << " != " << total_size
                << "total_size" << std::endl;
      exit(EXIT_FAILURE);
    }
    if constexpr (verbose)
    {
      std::cout << "# stddev\t" << gttl_stddev(size_vec.begin(),
                                               size_vec.end()) << std::endl;
    }
  }
  public:
  DNAEncodingMultiLength(const std::string &inputfilename)
      : total_size(0)
      , total_number_of_sequences(0)
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
      assert(sequence.size() < enc_vec.size());
      if (enc_vec[sequence.size()] == nullptr)
      {
        enc_vec[sequence.size()]
          = new ThisDNAEncodingForLength(sequence.size());
      }
      enc_vec[sequence.size()]->add(sequence);
    }
    size_t w_idx = 0;
    for (size_t idx = 0; idx < enc_vec.size(); idx++)
    {
      if (enc_vec[idx] != nullptr)
      {
        enc_vec[idx]->final_resize();
        assert(w_idx < idx);
        total_size += enc_vec[idx]->total_size_get();
        total_number_of_sequences += enc_vec[idx]->number_of_sequences_get();
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
  void prepare_view(size_t num_parts)
  {
    auto kv_vec = key_values_vector_get();
    if constexpr (verbose)
    {
      std::cout << "# original key values" << std::endl;
      key_values_show(kv_vec);
    }
    expanded_vec.clear();
    end_idx_of_part_vec.clear();
    if (num_parts == 1)
    {
      for (auto v : kv_vec)
      {
        expanded_vec.push_back(v);
      }
      end_idx_of_part_vec.push_back(expanded_vec.size());
      return;
    }
    const size_t mean = (total_size + num_parts - 1)/num_parts;
    if constexpr (verbose)
    {
      std::cout << "# num_parts\t" << num_parts << std::endl;
      std::cout << "# mean\t" << mean << std::endl;
    }
    for (auto &v : kv_vec)
    {
      const int size_factor = 15;
      const size_t s = std::get<1>(v) * std::get<2>(v);
      if (s < mean/size_factor)
      {
        expanded_vec.push_back(v);
      } else
      {
        const size_t mean_part_seq = std::max(size_t(10),
                                              mean/(size_factor *
                                                    std::get<2>(v)));
        size_t remain = std::get<1>(v);
        size_t sequence_number_offset = 0;
        while (remain >= mean_part_seq)
        {
          expanded_vec.push_back(std::make_tuple(std::get<0>(v),
                                                 mean_part_seq,
                                                 std::get<2>(v),
                                                 sequence_number_offset));
          remain -= mean_part_seq;
          sequence_number_offset += mean_part_seq;
        }
        if (remain > 0)
        {
          expanded_vec.push_back(std::make_tuple(std::get<0>(v),
                                                 remain,
                                                 std::get<2>(v),
                                                 sequence_number_offset));
        }
      }
    }
    if constexpr (verbose)
    {
      std::cout << "# expanded key values" << std::endl;
      key_values_show(expanded_vec);
    }
    divide_vector_evenly(num_parts);
  }

  class SplitViewIterator
  {
    const std::vector<ThisDNAEncodingForLength *> &enc_vec_ref;
    const KeyValuesType &expanded_vec_ref;
    const std::vector<size_t> &end_idx_of_part_vec_ref;
    const size_t end_of_part;
    size_t in_part_idx,
           current_enc_vec_idx;
    ThisDNAEncodingForLength *current_enc_vec_value;
    size_t num_units,
           sequence_number_offset;
    const uint64_t *current_units_ptr;
    size_t sequence_length,
           number_of_sequences;
    const uint64_t *units_end;
    bool exhausted;
    public:
    SplitViewIterator(const std::vector<ThisDNAEncodingForLength *>
                        &_enc_vec_ref,
                      const KeyValuesType &_expanded_vec_ref,
                      const std::vector<size_t> &_end_idx_of_part_vec_ref,
                      size_t part_idx,
                      bool _exhausted)
      : enc_vec_ref(_enc_vec_ref)
      , expanded_vec_ref(_expanded_vec_ref)
      , end_idx_of_part_vec_ref(_end_idx_of_part_vec_ref)
      , end_of_part(end_idx_of_part_vec_ref[part_idx])
      , in_part_idx(part_idx == 0 ? 0 : end_idx_of_part_vec_ref[part_idx - 1])
      , current_enc_vec_idx(std::get<0>(expanded_vec_ref[in_part_idx]))
      , current_enc_vec_value(enc_vec_ref[current_enc_vec_idx])
      , num_units(current_enc_vec_value->num_units_get())
      , sequence_number_offset(std::get<3>(expanded_vec_ref[in_part_idx]))
      , current_units_ptr(current_enc_vec_value->units_get() +
                          sequence_number_offset * num_units)
      , sequence_length(current_enc_vec_value->sequence_length_get())
      , number_of_sequences(std::get<1>(expanded_vec_ref[in_part_idx]))
      , units_end(current_units_ptr + number_of_sequences * num_units)
      , exhausted(_exhausted)
    {
      assert(part_idx < end_idx_of_part_vec_ref.size());
    }
    std::pair<const uint64_t *,size_t> operator *(void) const
    {
      return std::make_pair(current_units_ptr, sequence_length);
    }
    bool operator != (const SplitViewIterator& other) const noexcept
    {
      return exhausted != other.exhausted;
    }
    SplitViewIterator& operator++() /* prefix increment*/
    {
      assert(not exhausted);
      if (current_units_ptr + num_units < units_end)
      {
        current_units_ptr += num_units;
      } else
      {
        if (in_part_idx + 1 < end_of_part)
        {
          in_part_idx++;
          assert(in_part_idx < expanded_vec_ref.size());
          number_of_sequences = std::get<1>(expanded_vec_ref[in_part_idx]);
          current_enc_vec_idx = std::get<0>(expanded_vec_ref[in_part_idx]);
          assert(current_enc_vec_idx < enc_vec_ref.size());
          current_enc_vec_value = enc_vec_ref[current_enc_vec_idx];
          num_units = current_enc_vec_value->num_units_get();
          sequence_number_offset = std::get<3>(expanded_vec_ref[in_part_idx]);
          current_units_ptr = current_enc_vec_value->units_get()
                               + sequence_number_offset * num_units;
          sequence_length = current_enc_vec_value->sequence_length_get();
          units_end = current_units_ptr + number_of_sequences * num_units;
        } else
        {
          exhausted = true;
        }
      }
      return *this;
    }
  };

  auto begin(size_t part_idx = 0) const
  {
    if (expanded_vec.size() == 0 or end_idx_of_part_vec.size() == 0)
    {
      std::cerr << "need to call prepare_view before" << std::endl;
      exit(EXIT_FAILURE);
    }
    return SplitViewIterator(enc_vec,
                             expanded_vec,
                             end_idx_of_part_vec,
                             part_idx,
                             false);
  }
  auto end(void) const
  {
    return SplitViewIterator(enc_vec,
                             expanded_vec,
                             end_idx_of_part_vec,
                             0,
                             true);
  }
  size_t num_parts_get(void) const
  {
    return end_idx_of_part_vec.size();
  }
  size_t total_number_of_sequences_get(void) const
  {
    return total_number_of_sequences;
  }
  void verify_length_dist(std::map<size_t,size_t> &length_dist_map) const
  {
    size_t count_sum = 0;
    for (auto const& [len_as_key, count_as_value] : length_dist_map)
    {
      count_sum += count_as_value;
    }
    if (count_sum != total_number_of_sequences)
    {
      std::cerr << "count_sum = " << count_sum << " != "
                << total_number_of_sequences << " = total_number_of_sequences"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    for (auto &&dna_encoding : enc_vec)
    {
      const size_t this_length = dna_encoding->sequence_length_get();
      if (dna_encoding->number_of_sequences_get()
            != length_dist_map[this_length])
      {
        std::cerr << "dna_encoding->number_of_sequences_get() = "
                  << dna_encoding->number_of_sequences_get() << " != "
                  << length_dist_map[this_length] << " = "
                  << "length_dist_map[" << this_length << "]"
                  << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }
};
#endif
