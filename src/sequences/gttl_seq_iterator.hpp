#ifndef GTTL_SEQ_ITERATOR_HPP
#define GTTL_SEQ_ITERATOR_HPP
#include <iostream>
#include <string>
#include "utilities/gttl_line_iterator.hpp"

template<int buf_size>
class GttlSeqIterator
{
  struct Iterator
  {
    private:
      std::string &header,
                  &sequence;
      std::string *current_string = &header;
      GttlLineIterator<buf_size> &gttl_li;
      bool last_seq_was_processed,
           exhausted;
    public:
      Iterator(std::string &_header,
               std::string &_sequence,
               GttlLineIterator<buf_size> &_gttl_li,
               bool _exhausted) :
        header(_header),
        sequence(_sequence),
        gttl_li(_gttl_li),
        last_seq_was_processed(false),
        exhausted(_exhausted)
      {
      }
      std::pair<const std::string &,const std::string &> operator*()
      {
        if (!last_seq_was_processed)
        {
          header.clear();
          sequence.clear();
          current_string = &header;
          bool found_end = false;
          while (gttl_li.next(current_string))
          {
            if (current_string == &header)
            {
              current_string = &sequence;
            } else
            {
              sequence.pop_back();
              if (gttl_li.endofunit_get())
              {
                current_string = &header;
                found_end = true;
                break;
              }
            }
          }
          if (current_string == &sequence)
          {
            throw (std::ios_base::failure("last sequence has only header but "
                                          "no sequence content"));
          }
          if (!found_end)
          {
            exhausted = true;
          }
          last_seq_was_processed = true;
        }
        return {header,sequence};
      }
      Iterator& operator++() /* prefix increment*/
      {
        if (!gttl_li.more_lines())
        {
          exhausted = true;
        }
        last_seq_was_processed = false;
        return *this;
      }
      bool operator != (const Iterator& other) const
      {
        return exhausted != other.exhausted;
      }
  };
  private:
    std::string header{},
                sequence{};
    GttlLineIterator<buf_size> gttl_li;
  public:
    GttlSeqIterator(GttlFpType _in_fp) :
        gttl_li(GttlLineIterator<buf_size>(_in_fp))
    {
      gttl_li.separator_set('>');
    }
    Iterator begin()
    {
      return Iterator(header,sequence,gttl_li,false);
    }
    Iterator end()
    {
      return Iterator(header,sequence,gttl_li,true);
    }
};
#endif
