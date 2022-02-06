#ifndef GTTL_FASTQ_ITERATOR_HPP
#define GTTL_FASTQ_ITERATOR_HPP
#include <array>
#include <string>
#include <typeinfo>
#include <string_view>
#include "utilities/str_format.hpp"
#include "utilities/gttl_line_iterator.hpp"

struct FastQEntry
{
  std::array<std::string,4> seqbufs{};
  const std::string_view header_get(void) { return std::get<0>(seqbufs);}
  const std::string_view sequence_get(void) { return std::get<1>(seqbufs);}
  const std::string_view quality_get(void) { return std::get<3>(seqbufs);}
};

template<int buf_size>
class GttlFastQIterator
{
  struct Iterator
  {
    private:
      FastQEntry fastq_entry{};
      GttlLineIterator<buf_size> *gttl_li;
      bool last_seq_was_processed,
           input_exhausted;
    public:
      Iterator(GttlLineIterator<buf_size> *_gttl_li,
               bool _input_exhausted) :
        gttl_li(_gttl_li),
        last_seq_was_processed(false),
        input_exhausted(_input_exhausted) { }
      FastQEntry & operator*()
      {
        if (!last_seq_was_processed)
        {
          for (size_t idx = 0; idx < fastq_entry.seqbufs.size(); idx++)
          {
            fastq_entry.seqbufs[idx].clear();
          }
          int state = 0;
          bool found_end = false;
          while (gttl_li->next(&fastq_entry.seqbufs[state]))
          {
            fastq_entry.seqbufs[state].pop_back(); /* remove \n */
            if (state == 3)
            {
              found_end = true;
              break;
            }
            state++;
          }
          if (state != 3)
          {
            StrFormat msg(", line %lu: state=%d, corrupted sequence",
                          gttl_li->line_number_get()+1,state);
            throw msg.str(); /* check_err.py checked */
          }
          if (!found_end)
          {
            input_exhausted = true;
          }
          last_seq_was_processed = true;
        }
        if (fastq_entry.seqbufs[0].size() == 0)
        {
          StrFormat msg(", line %lu: XXcorrupted sequence",
                           gttl_li->line_number_get()+1);
          throw msg.str(); /* check_err.py checked */
        }
        return fastq_entry;
      }
      Iterator& operator++() /* prefix increment*/
      {
        if (!gttl_li->more_lines())
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
    GttlLineIterator<buf_size> *gttl_li;
    bool own_line_reader;
  public:
    GttlFastQIterator(GttlFpType _in_fp) :
      gttl_li(new GttlLineIterator<buf_size>(_in_fp)),
      own_line_reader(true)
    {
      gttl_li->separator_set('\n');
    }
    GttlFastQIterator(const char *inputfile) :
      gttl_li(new GttlLineIterator<buf_size>(inputfile)),
      own_line_reader(true)
    {
      gttl_li->separator_set('\n');
    }
    GttlFastQIterator(const std::string &inputfile) :
      gttl_li(new GttlLineIterator<buf_size>(inputfile.c_str())),
      own_line_reader(true)
    {
      gttl_li->separator_set('\n');
    }
    GttlFastQIterator(const std::vector<std::string> *inputfiles) :
      gttl_li(new GttlLineIterator<buf_size>(inputfiles)),
      own_line_reader(true)
    {
       gttl_li->separator_set('\n');
    }
    GttlFastQIterator(GttlLineIterator<buf_size> *_gttl_li) :
      gttl_li(_gttl_li),
      own_line_reader(false) {}
    ~GttlFastQIterator(void)
    {
      if (own_line_reader)
      {
        delete gttl_li;
      }
    }
    size_t line_number(void) const noexcept
    {
      return gttl_li->line_number_get();
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
