#ifndef MULTISEQ_FACTORY_HPP
#define MULTISEQ_FACTORY_HPP
#include <vector>
#include <string>
#include <cstdint>
#include <cmath>
#include <cstdio>
#include "sequences/gttl_multiseq.hpp"

class GttlMultiseqFactory
{
  private:
  std::vector<GttlMultiseq *> multiseq_vector;
  std::vector<size_t> seqnum_offset_vector;
  size_t num_parts,
         len_parts,
         num_sequences;;
  bool fastq_paired_input;
  public:
  GttlMultiseqFactory(const std::string &fastq_file0,
                      const std::string &fastq_file1,
                      size_t _num_sequences,
                      uint8_t padding_char,
                      bool short_header)
    : num_parts(0)
    , len_parts(0)
    , num_sequences(_num_sequences)
    , fastq_paired_input(true)
  {
    constexpr const int buf_size = 1 << 14;
    GttlLineIterator<buf_size> line_iterator0(fastq_file0.c_str()),
                               line_iterator1(fastq_file1.c_str());
    GttlFastQIterator<GttlLineIterator<buf_size>> fastq_it0(line_iterator0),
                                                  fastq_it1(line_iterator1);
    auto it0 = fastq_it0.begin();
    auto it1 = fastq_it1.begin();
    bool exhausted = false;

    while (!exhausted)
    {
      if (it0 == fastq_it0.end() || it1 == fastq_it1.end())
      {
        break;
      }
      GttlMultiseq *multiseq
        = new GttlMultiseq(true,padding_char); /* CONSTRUCTOR */

      for (size_t idx = 0; idx < num_sequences; idx++)
      {
        if (it0 == fastq_it0.end() || it1 == fastq_it1.end())
        {
          exhausted = true;
          break;
        }
        multiseq->append<true>((*it0).header_get(),(*it0).sequence_get(),
                               padding_char);
        multiseq->append<true>((*it1).header_get(),(*it1).sequence_get(),
                               padding_char);
        ++it0;
        ++it1;
      }
      if (multiseq->sequences_number_get() > 0)
      {
        if (short_header)
        {
          multiseq->short_header_cache_create<'|','|'>();
        }
        multiseq_vector.push_back(multiseq);
      } else
      {
        delete multiseq;
      }
    }
  }
  GttlMultiseqFactory(const std::string &inputfile,
                      size_t _num_parts,
                      size_t _len_parts,
                      size_t _num_sequences,
                      uint8_t padding_char,
                      bool short_header)
    : num_parts(_num_parts)
    , len_parts(_len_parts)
    , num_sequences(_num_sequences)
    , fastq_paired_input(false)
  {
    constexpr const int buf_size = 1 << 14;
    GttlSeqIterator<buf_size> fasta_it(inputfile.c_str());

    GttlMultiseq *multiseq
      = new GttlMultiseq(true,padding_char); /* CONSTRUCTOR */
    size_t seqnum = 0,
           current_part_number_of_units = 0,
           number_of_units_in_split = 0;
    if (num_parts > 0)
    {
      GttlSeqIterator<buf_size> fasta_it(inputfile.c_str());
      size_t sequences_total_length = 0;
      for (auto &&si : fasta_it)
      {
        sequences_total_length += si.sequence_get().size();
      }
      len_parts = sequences_total_length/num_parts;
    }
    if (len_parts > 0)
    {
      number_of_units_in_split = len_parts; /*total length of seq. in parts*/
    } else
    {
      number_of_units_in_split = num_sequences;
    }
    for (auto &&si : fasta_it)
    {
      if (current_part_number_of_units < number_of_units_in_split)
      {
        multiseq->append<true>(si.header_get(),si.sequence_get(),padding_char);
        if (len_parts > 0)
        {
          current_part_number_of_units += si.sequence_get().size();
        } else
        {
          assert(num_sequences > 0);
          current_part_number_of_units++;
        }
      } else
      {
        if (num_sequences > 0)
        {
          assert(current_part_number_of_units == number_of_units_in_split);
        }
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
          multiseq = new GttlMultiseq(true,padding_char); /* CONSTRUCTOR */
          multiseq->append<true>(si.header_get(),si.sequence_get(),
                                 padding_char);
        } else
        {
          delete multiseq;
          multiseq = nullptr;
        }
        current_part_number_of_units = 0;
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
  std::pair<GttlMultiseq *,size_t> at(size_t idx) const noexcept
  {
    assert(idx < multiseq_vector.size());
    size_t sequence_number_offset;
    if (len_parts > 0)
    {
      sequence_number_offset = idx == 0 ? 0 : seqnum_offset_vector[idx-1];
    } else
    {
      assert(num_sequences > 0);
      sequence_number_offset = idx * num_sequences;
    }
    return std::make_pair(multiseq_vector[idx],sequence_number_offset);
  }
  size_t size(void) const noexcept
  {
    return multiseq_vector.size();
  }
  bool is_fastq_paired_input(void) const noexcept
  {
    return fastq_paired_input;
  }
  void statistics(void) const noexcept
  {
    size_t total_length = 0, part_num = 0;
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
