#ifndef CHAR_RANGE_HPP
#define CHAR_RANGE_HPP
#include <cassert>
#include <cstddef>
#include <utility>

#define GTTL_CHAR_FINDER(VAR,INVERT,PTR)\
        if constexpr (INVERT)\
        {\
          if constexpr (forward)\
          {\
            VAR = char_finder.find_forward_not(PTR,end_ptr);\
          } else\
          {\
            VAR = char_finder.find_backward_not(PTR,end_ptr);\
          }\
        } else\
        {\
          if constexpr (forward)\
          {\
            VAR = char_finder.find_forward(PTR,end_ptr);\
          } else\
          {\
            VAR = char_finder.find_backward(PTR,end_ptr);\
          }\
        }

template<class CharFinder,const CharFinder &char_finder,bool forward,
         bool invert>
class GttlCharRange
{
  private:
  struct Iterator
  {
    private:
    const char *sequence, *curr_end, *end_ptr;
    size_t seqlen, range_start, range_length;
    bool exhausted;
    size_t ptr2difference(const char *a,const char *b)
    {
      if constexpr (forward)
      {
        assert(a>=b);
        return static_cast<size_t>(a - b);
      } else
      {
        assert(b>=a);
        return static_cast<size_t>(b - a);
      }
    }
    public:
      Iterator(const char *_sequence,size_t _seqlen,bool _exhausted) :
        curr_end(nullptr),
        seqlen(_seqlen),
        exhausted(_exhausted)
      {
        if (_sequence == nullptr)
        {
          sequence = end_ptr = nullptr;
        } else
        {
          if constexpr (forward)
          {
            sequence = _sequence;
            end_ptr = _sequence + _seqlen;
          } else
          {
            sequence = _sequence + _seqlen - 1;
            end_ptr = _sequence - 1;
          }
        }
        if (sequence != nullptr)
        {
          static constexpr const int step = forward ? 1 : -1;
          const char *curr_start;
          GTTL_CHAR_FINDER(curr_start,invert,sequence);
          if (curr_start != nullptr)
          {
            range_start = ptr2difference(curr_start,sequence);
            GTTL_CHAR_FINDER(curr_end,!invert,curr_start+step);
            if (curr_end == nullptr)
            {
              range_length = ptr2difference(end_ptr,curr_start);
              curr_start = nullptr;
            } else
            {
              range_length = ptr2difference(curr_end,curr_start);
            }
            if constexpr (!forward)
            {
              range_start = _seqlen - range_start - range_length;
            }
          } else
          {
            exhausted = true;
          }
        }
      }
      Iterator& operator++() /* prefix increment*/
      {
        if (curr_end == nullptr)
        {
          exhausted = true;
          return *this;
        }
        while (true)
        {
          static constexpr const int step = forward ? 1 : -1;
          assert(curr_end != nullptr);
          const char *curr_start;
          GTTL_CHAR_FINDER(curr_start,invert,curr_end+step);
          if (curr_start != nullptr)
          {
            range_start = ptr2difference(curr_start,sequence);
            GTTL_CHAR_FINDER(curr_end,!invert,curr_start+step);
            if (curr_end == nullptr)
            {
              range_length = ptr2difference(end_ptr,curr_start);
              if constexpr (!forward)
              {
                range_start = seqlen - range_start - range_length;
              }
              curr_start = nullptr;
              break;
            }
            range_length = ptr2difference(curr_end,curr_start);
            if constexpr (!forward)
            {
              range_start = seqlen - range_start - range_length;
            }
            break;
          } else
          {
            exhausted = true;
            break;
          }
        }
        return *this;
      }
      const std::pair<size_t,size_t> operator*(void) const
      {
        return {range_start,range_length};
      }
      bool operator != (const Iterator& other) const noexcept
      {
        return exhausted != other.exhausted;
      }
  };
  const char *sequence;
  size_t seqlen;
  public:
    GttlCharRange(const char *_sequence,size_t _seqlen) :
          sequence(_sequence),
          seqlen(_seqlen) {}
    [[nodiscard]] Iterator begin() const
    {
      return Iterator(sequence, seqlen, false);
    }
    [[nodiscard]] Iterator end() const { return Iterator(nullptr, 0, true); }
};

namespace char_range
{
  template<class CharFinder,const CharFinder &char_finder>
  using GttlForwardWildcardFinder
    = GttlCharRange<CharFinder,char_finder,true,true>;
  template<class CharFinder,const CharFinder &char_finder>
  using GttlBackwardWildcardFinder
    = GttlCharRange<CharFinder,char_finder,false,true>;
}
#endif
