#ifndef GTTL_FASTQ_ITERATOR_HPP
#define GTTL_FASTQ_ITERATOR_HPP
#include <array>
#include <string>
#include <typeinfo>
#include "utilities/str_format.hpp"
#include "utilities/gttl_line_iterator.hpp"

using FastQEntry = std::array<std::string,4>;

const std::string &fastq_header(const FastQEntry & fastq_entry)
{
  return std::get<0>(fastq_entry);
}

const std::string &fastq_sequence(const FastQEntry & fastq_entry)
{
  return std::get<1>(fastq_entry);
}

const std::string &fastq_quality(const FastQEntry & fastq_entry)
{
  return std::get<3>(fastq_entry);
}

template<int buf_size>
class GttlFastQIterator
{
  struct Iterator
  {
    private:
      FastQEntry seqbufs{};
      GttlLineIterator<buf_size> &gttl_li;
      bool last_seq_was_processed,
           input_exhausted;
    public:
      Iterator(GttlLineIterator<buf_size> &_gttl_li,
               bool _input_exhausted) :
        gttl_li(_gttl_li),
        last_seq_was_processed(false),
        input_exhausted(_input_exhausted) { }
      FastQEntry & operator*()
      {
        if (!last_seq_was_processed)
        {
          for (size_t idx = 0; idx < seqbufs.size(); idx++)
          {
            seqbufs[idx].clear();
          }
          int state = 0;
          bool found_end = false;
          while (gttl_li.next(&seqbufs[state]))
          {
            seqbufs[state].pop_back();
            if (state == 3)
            {
              found_end = true;
              break;
            }
            state++;
          }
          if (state != 3)
          {
            StrFormat msg(", line %lu: state=%d,YYcorrupted sequence",
                           gttl_li.line_number_get()+1,state);
            throw msg.str(); /* check_err.py checked */
          }
          if (!found_end)
          {
            input_exhausted = true;
          }
          last_seq_was_processed = true;
        }
        if (seqbufs[0].size() == 0)
        {
          StrFormat msg(", line %lu: XXcorrupted sequence",
                           gttl_li.line_number_get()+1);
          throw msg.str(); /* check_err.py checked */
        }
        return seqbufs;
      }
      Iterator& operator++() /* prefix increment*/
      {
        if (!gttl_li.more_lines())
        {
          input_exhausted = true;
        }
        last_seq_was_processed = false;
        return *this;
      }
      bool operator != (const Iterator& other) const
      {
        return input_exhausted != other.input_exhausted;
      }
      bool operator == (const Iterator& other) const
      {
        return input_exhausted == other.input_exhausted;
      }
  };
  private:
    GttlLineIterator<buf_size> gttl_li;
  public:
    GttlFastQIterator(GttlFpType _in_fp) :
        gttl_li(GttlLineIterator<buf_size>(_in_fp))
    {
      gttl_li.separator_set('\n');
    }
    GttlFastQIterator(const char *inputfile) :
        gttl_li(GttlLineIterator<buf_size>(inputfile))
    {
      gttl_li.separator_set('\n');
    }
    GttlFastQIterator(const std::string &inputfile) :
        gttl_li(GttlLineIterator<buf_size>(inputfile.c_str()))
    {
      gttl_li.separator_set('\n');
    }
    GttlFastQIterator(const std::vector<std::string> *inputfiles) :
      gttl_li(GttlLineIterator<buf_size>(inputfiles))
    {
       gttl_li.separator_set('\n');
    }
    size_t line_number(void) const noexcept
    {
      return gttl_li.line_number_get();
    }
    Iterator begin()
    {
      return Iterator(gttl_li,false);
    }
    Iterator end()
    {
      return Iterator(gttl_li,true);
    }
};

#endif
