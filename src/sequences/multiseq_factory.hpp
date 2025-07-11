#ifndef MULTISEQ_FACTORY_HPP
#define MULTISEQ_FACTORY_HPP
#include <cassert>
#include <utility>
#include <vector>
#include <string>
#include <cstdint>
#include <cmath>
#include <cstdio>
#include "sequences/gttl_fasta_generator.hpp"
#include "sequences/gttl_fastq_generator.hpp"
#include "sequences/gttl_multiseq.hpp"

class GttlMultiseqFactory
{
  private:
  static constexpr const bool store_header = true;
  static constexpr const bool store_sequence = true;
  static constexpr const int buf_size = 1 << 14;
  std::vector<GttlMultiseq *> multiseq_vector;
  std::vector<size_t> seqnum_offset_vector;
  size_t num_sequences;
  bool fastq_paired_input;
  [[nodiscard]] size_t
  fastq_file_total_length_get(const std::string &inputfile) const
  {
    GttlFastQGenerator<buf_size> fastq_it(inputfile.c_str());
    size_t sequences_total_length = 0;
    for (const auto *it : fastq_it)
    {
      sequences_total_length += it->sequence_get().size();
    }
    return sequences_total_length;
  }
  public:
  GttlMultiseqFactory(const std::string &fastq_file0,
                      const std::string &fastq_file1,
                      size_t num_parts,
                      size_t len_parts,
                      size_t _num_sequences,
                      uint8_t padding_char,
                      bool short_header)
    : num_sequences(_num_sequences)
    , fastq_paired_input(true)
  {
    GttlMultiseq *multiseq
      = new GttlMultiseq(store_sequence, padding_char);/*CONSTRUCTOR */
    if (num_parts > 0)
    {
      const size_t sequences_total_length
        = fastq_file_total_length_get(fastq_file0) +
          fastq_file_total_length_get(fastq_file1);
      assert(len_parts == 0);
      len_parts = sequences_total_length/num_parts;
    }
    size_t number_of_units_in_split;
    if (len_parts > 0)
    {
      number_of_units_in_split = len_parts; /*total length of seq. in parts*/
    } else
    {
      number_of_units_in_split = num_sequences;
    }
    size_t seqnum = 0;
    size_t current_part_number_of_units = 0;
    GttlFastQGenerator<buf_size> fastq_it0(fastq_file0.c_str());
    GttlFastQGenerator<buf_size> fastq_it1(fastq_file1.c_str());
    auto it0 = fastq_it0.begin();
    auto it1 = fastq_it1.begin();
    while(it0 != fastq_it0.end() and it1 != fastq_it1.end())
    {
      if (current_part_number_of_units < number_of_units_in_split)
      {
        auto sequence0 = (*it0)->sequence_get();
        multiseq->append((*it0)->header_get(),sequence0,
                         store_header, store_sequence,padding_char);
        auto sequence1 = (*it1)->sequence_get();
        multiseq->append((*it1)->header_get(),sequence1,
                         store_header, store_sequence, padding_char);
        current_part_number_of_units
          += (len_parts > 0 ? (sequence0.size() + sequence1.size())
                            : 2);
      } else
      {
        assert(len_parts > 0 or
               current_part_number_of_units == number_of_units_in_split);
        if (multiseq->sequences_number_get() > 0)
        {
          if (short_header)
          {
            multiseq->short_header_cache_create<'|','|'>();
          }
          multiseq_vector.push_back(multiseq);
          if (len_parts > 0)
          {
            seqnum_offset_vector.push_back(seqnum);
          }
          multiseq = new GttlMultiseq(store_sequence,
                                      padding_char); /* CONSTRUCTOR */
          auto sequence0 = (*it0)->sequence_get();
          multiseq->append((*it0)->header_get(),sequence0,
                           store_header,store_sequence,padding_char);
          auto sequence1 = (*it1)->sequence_get();
          multiseq->append((*it1)->header_get(),sequence1,
                           store_sequence,store_header,padding_char);
          current_part_number_of_units
            = (len_parts > 0 ? (sequence0.size() + sequence1.size())
                             : 2);
        } else
        {
          delete multiseq;
          multiseq = nullptr;
          current_part_number_of_units = 0;
        }
      }
      ++it0;
      ++it1;
      seqnum += 2;
    }
    if (multiseq != nullptr and multiseq->sequences_number_get() > 0)
    {
      if (short_header)
      {
        multiseq->short_header_cache_create<'|','|'>();
      }
      multiseq_vector.push_back(multiseq);
      if (len_parts > 0)
      {
        seqnum_offset_vector.push_back(seqnum);
      }
    }
  }
  GttlMultiseqFactory(const std::string &inputfile,
                      size_t num_parts,
                      size_t len_parts,
                      size_t _num_sequences,
                      uint8_t padding_char,
                      bool short_header)
    : num_sequences(_num_sequences)
    , fastq_paired_input(false)
  {
    constexpr const int buf_size = 1 << 14;
    GttlFastAGenerator<buf_size> fasta_it(inputfile.c_str());

    GttlMultiseq *multiseq
      = new GttlMultiseq(store_sequence, padding_char); /* CONSTRUCTOR */
    if (num_parts > 0)
    {
      GttlFastAGenerator<buf_size> fasta_it(inputfile.c_str());
      size_t sequences_total_length = 0;
      for (auto &&si : fasta_it)
      {
        sequences_total_length += si->sequence_get().size();
      }
      assert(len_parts == 0);
      len_parts = sequences_total_length/num_parts;
    }
    size_t number_of_units_in_split;
    if (len_parts > 0)
    {
      number_of_units_in_split = len_parts; /*total length of seq. in parts*/
    } else
    {
      number_of_units_in_split = num_sequences;
    }
    size_t seqnum = 0;
    size_t current_part_number_of_units = 0;
    for (auto &&si : fasta_it)
    {
      if (current_part_number_of_units < number_of_units_in_split)
      {
        multiseq->append(si->header_get(),
                         si->sequence_get(),
                         store_header,store_sequence,padding_char);
        current_part_number_of_units
          += (len_parts > 0 ? si->sequence_get().size() : 1);
      } else
      {
        assert(len_parts > 0 or
               current_part_number_of_units == number_of_units_in_split);
        if (multiseq->sequences_number_get() > 0)
        {
          if (short_header)
          {
            multiseq->short_header_cache_create<'|','|'>();
          }
          multiseq_vector.push_back(multiseq);
          if (len_parts > 0)
          {
            seqnum_offset_vector.push_back(seqnum);
          }
          multiseq = new GttlMultiseq(store_sequence,
                                      padding_char); /* CONSTRUCTOR */
          multiseq->append(si->header_get(),si->sequence_get(),
                           store_header,store_sequence,padding_char);
          current_part_number_of_units
            = (len_parts > 0 ? si->sequence_get().size() : 1);
        } else
        {
          delete multiseq;
          multiseq = nullptr;
          current_part_number_of_units = 0;
        }
      }
      seqnum++;
    }
    if (multiseq != nullptr and multiseq->sequences_number_get() > 0)
    {
      if (short_header)
      {
        multiseq->short_header_cache_create<'|','|'>();
      }
      multiseq_vector.push_back(multiseq);
      if (len_parts > 0)
      {
        seqnum_offset_vector.push_back(seqnum);
      }
    }
  }
  ~GttlMultiseqFactory(void)
  {
    for (auto &&ms : multiseq_vector)
    {
      delete ms;
    }
  }
  [[nodiscard]] std::pair<GttlMultiseq *, size_t> at(size_t idx) const noexcept
  {
    assert(idx < multiseq_vector.size());
    size_t sequence_number_offset;
    if (seqnum_offset_vector.size() > 0)
    {
      sequence_number_offset = idx == 0 ? 0 : seqnum_offset_vector[idx-1];
    } else
    {
      assert(num_sequences > 0);
      sequence_number_offset = idx * num_sequences;
    }
    return std::make_pair(multiseq_vector[idx],sequence_number_offset);
  }
  [[nodiscard]] size_t size(void) const noexcept
  {
    return multiseq_vector.size();
  }
  [[nodiscard]] bool is_fastq_paired_input(void) const noexcept
  {
    return fastq_paired_input;
  }
  void statistics(void) const noexcept
  {
    size_t total_length = 0;
    size_t part_num = 0;
    for (auto &&ms_ptr : multiseq_vector)
    {
      total_length += ms_ptr->sequences_total_length_get();
      printf("%zu\t%zu\n",part_num,ms_ptr->sequences_total_length_get());
      part_num++;
    }
    const size_t mean = total_length/this->size();
    double squared_difference = 0;
    printf("# total_length\t%zu\n",total_length);
    printf("# mean length\t%zu\n",mean);
    for (auto &&ms_ptr : multiseq_vector)
    {
      const double diff
        = static_cast<double>(ms_ptr->sequences_total_length_get()) - mean;
      squared_difference += (diff * diff);
    }
    const double variance = squared_difference / this->size();
    printf("# stddev of size of parts\t%.2e\n",std::sqrt(variance));
  }
  void sequence_output(size_t sequence_output_width) const noexcept
  {
    for (auto &&ms_ptr : multiseq_vector)
    {
      ms_ptr->show(sequence_output_width, false);
    }
  }
};
#endif
