#ifndef MULTISEQ_FACTORY_HPP
#define MULTISEQ_FACTORY_HPP
#include <vector>
#include <string>
#include <cstdint>
#include "sequences/gttl_multiseq.hpp"

class GttlMultiseqFactory
{
  private:
  std::vector<GttlMultiseq *> multiseq_vector;
  size_t number_of_sequences_in_split;
  bool fastq_paired_input;
  public:
  GttlMultiseqFactory(const std::string &fastq_file0,
                      const std::string &fastq_file1,
                      size_t _number_of_sequences_in_split,
                      uint8_t padding_char,
                      bool short_header)
    : number_of_sequences_in_split(_number_of_sequences_in_split)
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

      for (size_t idx = 0; idx < number_of_sequences_in_split; idx++)
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
                      size_t _number_of_sequences_in_split,
                      uint8_t padding_char,
                      bool short_header)
    : number_of_sequences_in_split(_number_of_sequences_in_split)
    , fastq_paired_input(false)
  {
    constexpr const int buf_size = 1 << 14;
    GttlSeqIterator<buf_size> fasta_it(inputfile.c_str());

    size_t current_part_number_of_sequences = 0;
    GttlMultiseq *multiseq
      = new GttlMultiseq(true,padding_char); /* CONSTRUCTOR */
    for (auto &&si : fasta_it)
    {
      if (current_part_number_of_sequences < number_of_sequences_in_split)
      {
        multiseq->append<true>(si.header_get(),si.sequence_get(),padding_char);
        current_part_number_of_sequences++;
      } else
      {
        assert(current_part_number_of_sequences
                 == number_of_sequences_in_split);
        if (multiseq->sequences_number_get() > 0)
        {
          if (short_header)
          {
            multiseq->short_header_cache_create<'|','|'>();
          }
          multiseq_vector.push_back(multiseq);
          multiseq = new GttlMultiseq(true,padding_char); /* CONSTRUCTOR */
        } else
        {
          delete multiseq;
          multiseq = nullptr;
        }
        current_part_number_of_sequences = 0;
      }
    }
    if (multiseq != nullptr and multiseq->sequences_number_get() > 0)
    {
      if (short_header)
      {
        multiseq->short_header_cache_create<'|','|'>();
      }
      multiseq_vector.push_back(multiseq);
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
    const size_t sequence_number_offset = idx * number_of_sequences_in_split;
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
};
#endif
