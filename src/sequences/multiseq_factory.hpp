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
  size_t split_size;
  public:
  GttlMultiseqFactory(const std::string &filename0,
                      const std::string &filename1,
                      size_t _split_size,
                      uint8_t padding_char,
                      bool short_header)
    : multiseq_vector({})
    , split_size(_split_size)
  {
    constexpr const int buf_size = 1 << 14;
    GttlLineIterator<buf_size> line_iterator0(filename0.c_str()),
                               line_iterator1(filename1.c_str());
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

      for (size_t idx = 0; idx < split_size; idx++)
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
    return std::make_pair(multiseq_vector[idx],idx * split_size);
  }
  size_t size(void) const noexcept
  {
    return multiseq_vector.size();
  }
};
#endif
