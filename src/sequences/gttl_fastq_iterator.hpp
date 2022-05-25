#ifndef GTTL_FASTQ_ITERATOR_HPP
#define GTTL_FASTQ_ITERATOR_HPP
#include <array>
#include <string>
#include <typeinfo>
#include <type_traits>
#include <string_view>
#include "utilities/str_format.hpp"
#include "utilities/gttl_line_iterator.hpp"

template<int buf_size>
class GttlFastQIterator
{
  using LineIterator = GttlLineIterator<buf_size>;
  template<class StringStoreType>
  struct Iterator
  {
    struct FastQEntry
    {
      std::array<StringStoreType,4> seqbufs{};
      private:
      template<int state>
      auto access(void)
      {
        if constexpr (std::is_same_v<StringStoreType, std::string>)
        {
          return static_cast<const std::string_view>(std::get<state>(seqbufs));
        } else
        {
          static_assert(std::is_same_v<StringStoreType,ModifiableStringView>);
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
      LineIterator *gttl_li;
      bool last_seq_was_processed,
           input_exhausted;
    public:
      Iterator(LineIterator *_gttl_li,
               bool _input_exhausted)
        : gttl_li(_gttl_li)
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
          StrFormat msg(", line %lu: corrupted sequence",
                           gttl_li->line_number_get()+1);
          throw msg.str(); /* check_err.py checked */
        }
        return fastq_entry;
      }
      Iterator<StringStoreType>& operator++() /* prefix increment*/
      {
        if (!gttl_li->more_lines())
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
    LineIterator *gttl_li;
    bool own_line_iterator;
  public:
    GttlFastQIterator(GttlFpType _in_fp)
      : gttl_li(new LineIterator(_in_fp))
      , own_line_iterator(true)
    {
      gttl_li->separator_set('\n');
    }
    GttlFastQIterator(const char *inputfile)
      : gttl_li(new LineIterator(inputfile))
      , own_line_iterator(true)
    {
      gttl_li->separator_set('\n');
    }
    GttlFastQIterator(const std::string &inputfile)
      : gttl_li(new LineIterator(inputfile.c_str()))
      , own_line_iterator(true)
    {
      gttl_li->separator_set('\n');
    }
    GttlFastQIterator(const std::vector<std::string> *inputfiles)
      : gttl_li(new LineIterator(inputfiles))
      , own_line_iterator(true)
    {
      gttl_li->separator_set('\n');
    }
    GttlFastQIterator(LineIterator *_gttl_li)
      : gttl_li(_gttl_li)
      , own_line_iterator(false)
      {}
    GttlFastQIterator(const char *input_string,size_t len)
      : gttl_li(new LineIterator(input_string,len))
      , own_line_iterator(false)
    {
      gttl_li->separator_set('\n');
    }
    ~GttlFastQIterator(void)
    {
      if (own_line_iterator)
      {
        delete gttl_li;
      }
    }
    size_t line_number(void) const noexcept
    {
      return gttl_li->line_number_get();
    }
    auto begin(void)
    {
      if constexpr (buf_size == 0)
      {
        return Iterator<ModifiableStringView>(gttl_li,false);
      } else
      {
        return Iterator<std::string>(gttl_li,false);
      }
    }
    auto end(void)
    {
      if constexpr (buf_size == 0)
      {
        return Iterator<ModifiableStringView>(gttl_li,true);
      } else
      {
        return Iterator<std::string>(gttl_li,true);
      }
    }
};
#endif
