#ifndef GTTL_FASTQ_ITERATOR_HPP
#define GTTL_FASTQ_ITERATOR_HPP
#include <array>
#include <cstddef>
#include <ios>
#include <string>
#include <type_traits>
#include <string_view>
#include "utilities/str_format.hpp"
#include "utilities/gttl_line_iterator.hpp"

template<class LineIterator>
class GttlFastQIterator
{
  template<class StringStoreType>
  struct Iterator
  {
    struct FastQEntry
    {
      std::array<StringStoreType,4> seqbufs{};
      private:
      template<int state>
      auto access(void) const noexcept
      {
        if constexpr (std::is_same_v<StringStoreType, std::string>)
        {
          return static_cast<const std::string_view>(std::get<state>(seqbufs));
        } else
        {
          return std::string_view(std::get<state>(seqbufs).data(),
                                  std::get<state>(seqbufs).size());
        }
      }
      public:
      auto header_get(void) {return access<0>(); }
      auto sequence_get(void) {return access<1>(); }
      auto quality_get(void) {return access<3>(); }
    };
    private:
      FastQEntry fastq_entry{};
      LineIterator &line_iterator;
      bool last_seq_was_processed,
           input_exhausted;
    public:
      Iterator(LineIterator &_line_iterator,
               bool _input_exhausted)
        : line_iterator(_line_iterator)
        , last_seq_was_processed(false)
        , input_exhausted(_input_exhausted)
        {}
      FastQEntry &operator*()
      {
        if (!last_seq_was_processed)
        {
          for (size_t idx = 0; idx < fastq_entry.seqbufs.size(); idx++)
          {
            fastq_entry.seqbufs[idx].clear();
          }
          int state = 0;
          bool found_end = false;
          while (line_iterator.next(&fastq_entry.seqbufs[state]))
          {
            fastq_entry.seqbufs[state].pop_back(); /* remove trailing \n */
            if (state == 3)
            {
              found_end = true;
              break;
            }
            state++;
          }
          if (state != 3)
          {
            StrFormat msg(", line %zu: state=%d, corrupted sequence",
                          line_iterator.line_number_get()+1,state);
            throw std::ios_base::failure(msg.str()); /* check_err.py checked */
          }
          if (!found_end)
          {
            input_exhausted = true;
          }
          last_seq_was_processed = true;
        }
        if (fastq_entry.seqbufs[0].size() == 0)
        {
          StrFormat msg(", line %zu: corrupted sequence",
                           line_iterator.line_number_get()+1);
          throw std::ios_base::failure(msg.str()); /* check_err.py checked */
        }
        return fastq_entry;
      }
      Iterator<StringStoreType>& operator++() /* prefix increment*/
      {
        if (!line_iterator.more_lines())
        {
          input_exhausted = true;
        }
        last_seq_was_processed = false;
        return *this;
      }
      bool operator != (const Iterator<StringStoreType>& other) const
      {
        return input_exhausted != other.input_exhausted;
      }
      bool operator == (const Iterator<StringStoreType>& other) const
      {
        return input_exhausted == other.input_exhausted;
      }
  };
  private:
    LineIterator &line_iterator;
  public:
    static constexpr bool is_fastq_iterator = true;
    GttlFastQIterator(LineIterator &_line_iterator)
      : line_iterator(_line_iterator)
      {}
    void reset(void)
    {
      line_iterator.reset();
    }
    size_t line_number(void) const noexcept
    {
      return line_iterator->line_number_get();
    }
    auto begin(void)
    {
      if constexpr (LineIterator::this_buf_size > 0)
      {
        return Iterator<std::string>(line_iterator,false);
      } else
      {
        return Iterator<LineIteratorSubstring>(line_iterator,false);
      }
    }
    auto end(void)
    {
      if constexpr (LineIterator::this_buf_size > 0)
      {
        return Iterator<std::string>(line_iterator,true);
      } else
      {
        return Iterator<LineIteratorSubstring>(line_iterator,true);
      }
    }
};
#endif
