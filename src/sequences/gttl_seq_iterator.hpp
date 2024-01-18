#ifndef GTTL_SEQ_ITERATOR_HPP
#define GTTL_SEQ_ITERATOR_HPP
#include <iostream>
#include <string>
#include <tuple>
#include "utilities/gttl_line_iterator.hpp"

struct SequenceEntry
{
  std::string header{},
              sequence{};
  SequenceEntry(void) {}
  SequenceEntry(std::string _header,std::string _sequence)
    : header(_header)
    , sequence(_sequence)
  {}
  const std::string_view header_get(void) { return header;}
  const std::string_view sequence_get(void) {return sequence;}
};

template<int buf_size>
class GttlSeqIterator
{
  struct Iterator
  {
    private:
      SequenceEntry sequence_entry;
      std::string *current_string;
      GttlLineIterator<buf_size> &gttl_li;
      bool last_seq_was_processed,
           input_exhausted;
    public:
      Iterator(GttlLineIterator<buf_size> &_gttl_li,
               bool _input_exhausted)
        : sequence_entry({})
        , current_string(&sequence_entry.header)
        , gttl_li(_gttl_li)
        , last_seq_was_processed(false)
        , input_exhausted(_input_exhausted)
      {
      }
      SequenceEntry &operator*()
      {
        if (!last_seq_was_processed)
        {
          sequence_entry.header.clear();
          sequence_entry.sequence.clear();
          current_string = &sequence_entry.header;
          bool found_end = false;
          while (gttl_li.next(current_string))
          {
            current_string->pop_back(); /* remove \n */
            if (current_string == &sequence_entry.header)
            {
              current_string = &sequence_entry.sequence;
            } else
            {
              if (gttl_li.endofunit_get())
              {
                current_string = &sequence_entry.header;
                found_end = true;
                break;
              }
            }
          }
          if (current_string == &sequence_entry.sequence)
          {
            StrFormat msg(", line %lu: corrupted sequence",
                           gttl_li.line_number_get()+1);
            throw msg.str(); /* check_err.py checked */
          }
          if (!found_end)
          {
            input_exhausted = true;
          }
          last_seq_was_processed = true;
        }
        if (sequence_entry.sequence.size() == 0 ||
            sequence_entry.sequence[0] == '>')
        {
          StrFormat msg(", line %lu: corrupted sequence",
                           gttl_li.line_number_get()+1);
          throw msg.str(); /* check_err.py checked */
        }
        return sequence_entry;
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
  };
  private:
    GttlLineIterator<buf_size> gttl_li;
  public:
    GttlSeqIterator(GttlFpType _in_fp) :
        gttl_li(GttlLineIterator<buf_size>(_in_fp))
    {
      gttl_li.separator_set('>');
    }
    GttlSeqIterator(const char *inputfile) :
        gttl_li(GttlLineIterator<buf_size>(inputfile))
    {
      gttl_li.separator_set('>');
    }
    GttlSeqIterator(const std::vector<std::string> *inputfiles) :
      gttl_li(GttlLineIterator<buf_size>(inputfiles))
    {
      gttl_li.separator_set('>');
    }
    GttlSeqIterator(const std::string_view &multi_fasta_part)
      : gttl_li(GttlLineIterator<0>(multi_fasta_part.data(),
                                    multi_fasta_part.size()))
    {
      gttl_li.separator_set('>');
    }
    size_t line_number(void) const noexcept
    {
      return gttl_li.line_number_get();
    }
    Iterator begin()
    {
      gttl_li.reset();
      return Iterator(gttl_li,false);
    }
    Iterator end()
    {
      return Iterator(gttl_li,true);
    }
};
#endif
