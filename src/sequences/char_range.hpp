#ifndef CHAR_RANGE_HPP
#define CHAR_RANGE_HPP
#include <array>
#include <iostream>
#include <cstdbool>
#include <cstdint>

template <const char *charset>
static constexpr bool is_member_in_charset_rec(int i,char cc)
{
  return charset[i] == '\0'
           ? false
           : (charset[i] == cc ? true
                               : is_member_in_charset_rec<charset>(i+1,cc));
}

template <const char *charset,int cc>
static constexpr bool is_member_in_charset(void)
{
  return is_member_in_charset_rec<charset>(0,static_cast<char>(cc));
}

template<const char *charset>
class CharFinder
{
  static constexpr const bool in_charset[] =
  {
    /* 0 */ false,
    /* 1 */ false,
    /* 2 */ false,
    /* 3 */ false,
    /* 4 */ false,
    /* 5 */ false,
    /* 6 */ false,
    /* 7 */ false,
    /* 8 */ false,
    /* 9 */ false,
    /* 10 */ false,
    /* 11 */ false,
    /* 12 */ false,
    /* 13 */ false,
    /* 14 */ false,
    /* 15 */ false,
    /* 16 */ false,
    /* 17 */ false,
    /* 18 */ false,
    /* 19 */ false,
    /* 20 */ false,
    /* 21 */ false,
    /* 22 */ false,
    /* 23 */ false,
    /* 24 */ false,
    /* 25 */ false,
    /* 26 */ false,
    /* 27 */ false,
    /* 28 */ false,
    /* 29 */ false,
    /* 30 */ false,
    /* 31 */ false,
    /* 32 */ false,
    /* 33 */ is_member_in_charset<charset,33>(),
    /* 34 */ is_member_in_charset<charset,34>(),
    /* 35 */ is_member_in_charset<charset,35>(),
    /* 36 */ is_member_in_charset<charset,36>(),
    /* 37 */ is_member_in_charset<charset,37>(),
    /* 38 */ is_member_in_charset<charset,38>(),
    /* 39 */ is_member_in_charset<charset,39>(),
    /* 40 */ is_member_in_charset<charset,40>(),
    /* 41 */ is_member_in_charset<charset,41>(),
    /* 42 */ is_member_in_charset<charset,42>(),
    /* 43 */ is_member_in_charset<charset,43>(),
    /* 44 */ is_member_in_charset<charset,44>(),
    /* 45 */ is_member_in_charset<charset,45>(),
    /* 46 */ is_member_in_charset<charset,46>(),
    /* 47 */ is_member_in_charset<charset,47>(),
    /* 48 */ is_member_in_charset<charset,48>(),
    /* 49 */ is_member_in_charset<charset,49>(),
    /* 50 */ is_member_in_charset<charset,50>(),
    /* 51 */ is_member_in_charset<charset,51>(),
    /* 52 */ is_member_in_charset<charset,52>(),
    /* 53 */ is_member_in_charset<charset,53>(),
    /* 54 */ is_member_in_charset<charset,54>(),
    /* 55 */ is_member_in_charset<charset,55>(),
    /* 56 */ is_member_in_charset<charset,56>(),
    /* 57 */ is_member_in_charset<charset,57>(),
    /* 58 */ is_member_in_charset<charset,58>(),
    /* 59 */ is_member_in_charset<charset,59>(),
    /* 60 */ is_member_in_charset<charset,60>(),
    /* 61 */ is_member_in_charset<charset,61>(),
    /* 62 */ is_member_in_charset<charset,62>(),
    /* 63 */ is_member_in_charset<charset,63>(),
    /* 64 */ is_member_in_charset<charset,64>(),
    /* 65 */ is_member_in_charset<charset,65>(),
    /* 66 */ is_member_in_charset<charset,66>(),
    /* 67 */ is_member_in_charset<charset,67>(),
    /* 68 */ is_member_in_charset<charset,68>(),
    /* 69 */ is_member_in_charset<charset,69>(),
    /* 70 */ is_member_in_charset<charset,70>(),
    /* 71 */ is_member_in_charset<charset,71>(),
    /* 72 */ is_member_in_charset<charset,72>(),
    /* 73 */ is_member_in_charset<charset,73>(),
    /* 74 */ is_member_in_charset<charset,74>(),
    /* 75 */ is_member_in_charset<charset,75>(),
    /* 76 */ is_member_in_charset<charset,76>(),
    /* 77 */ is_member_in_charset<charset,77>(),
    /* 78 */ is_member_in_charset<charset,78>(),
    /* 79 */ is_member_in_charset<charset,79>(),
    /* 80 */ is_member_in_charset<charset,80>(),
    /* 81 */ is_member_in_charset<charset,81>(),
    /* 82 */ is_member_in_charset<charset,82>(),
    /* 83 */ is_member_in_charset<charset,83>(),
    /* 84 */ is_member_in_charset<charset,84>(),
    /* 85 */ is_member_in_charset<charset,85>(),
    /* 86 */ is_member_in_charset<charset,86>(),
    /* 87 */ is_member_in_charset<charset,87>(),
    /* 88 */ is_member_in_charset<charset,88>(),
    /* 89 */ is_member_in_charset<charset,89>(),
    /* 90 */ is_member_in_charset<charset,90>(),
    /* 91 */ is_member_in_charset<charset,91>(),
    /* 92 */ is_member_in_charset<charset,92>(),
    /* 93 */ is_member_in_charset<charset,93>(),
    /* 94 */ is_member_in_charset<charset,94>(),
    /* 95 */ is_member_in_charset<charset,95>(),
    /* 96 */ is_member_in_charset<charset,96>(),
    /* 97 */ is_member_in_charset<charset,97>(),
    /* 98 */ is_member_in_charset<charset,98>(),
    /* 99 */ is_member_in_charset<charset,99>(),
    /* 100 */ is_member_in_charset<charset,100>(),
    /* 101 */ is_member_in_charset<charset,101>(),
    /* 102 */ is_member_in_charset<charset,102>(),
    /* 103 */ is_member_in_charset<charset,103>(),
    /* 104 */ is_member_in_charset<charset,104>(),
    /* 105 */ is_member_in_charset<charset,105>(),
    /* 106 */ is_member_in_charset<charset,106>(),
    /* 107 */ is_member_in_charset<charset,107>(),
    /* 108 */ is_member_in_charset<charset,108>(),
    /* 109 */ is_member_in_charset<charset,109>(),
    /* 110 */ is_member_in_charset<charset,110>(),
    /* 111 */ is_member_in_charset<charset,111>(),
    /* 112 */ is_member_in_charset<charset,112>(),
    /* 113 */ is_member_in_charset<charset,113>(),
    /* 114 */ is_member_in_charset<charset,114>(),
    /* 115 */ is_member_in_charset<charset,115>(),
    /* 116 */ is_member_in_charset<charset,116>(),
    /* 117 */ is_member_in_charset<charset,117>(),
    /* 118 */ is_member_in_charset<charset,118>(),
    /* 119 */ is_member_in_charset<charset,119>(),
    /* 120 */ is_member_in_charset<charset,120>(),
    /* 121 */ is_member_in_charset<charset,121>(),
    /* 122 */ is_member_in_charset<charset,122>(),
    /* 123 */ is_member_in_charset<charset,123>(),
    /* 124 */ false,
    /* 125 */ is_member_in_charset<charset,125>(),
    /* 126 */ is_member_in_charset<charset,126>(),
    /* 127 */ false
  };
  private:
  template<bool ref_value>
  const char *find_generic(const char *s,const char *endptr) const
  {
    for (const char *sptr = s; sptr < endptr; ++sptr)
    {
      if (in_charset[static_cast<int>(*sptr)] == ref_value)
      {
        return sptr;
      }
    }
    return nullptr;
  }
  public:
  bool is_member(char cc) const
  {
    return in_charset[static_cast<int>(cc)];
  }
  const char *find(const char *s,const char *endptr) const
  {
    return find_generic<true>(s,endptr);
  }
  const char *find_not(const char *s,const char *endptr) const
  {
    return find_generic<false>(s,endptr);
  }
};

template<const char *charset>
class GttlCharRange
{
  private:
  struct Iterator
  {
    private:
    static constexpr const CharFinder<charset> char_finder{};
    const char *sequence, *curr_end, *end_ptr;
    size_t range_start, range_length;
    bool exhausted;
    public:
      Iterator(const char *_sequence,size_t seqlen,bool _exhausted) :
        sequence(_sequence),
        curr_end(nullptr),
        end_ptr(sequence == nullptr ? nullptr : (sequence + seqlen)),
        exhausted(_exhausted)
      {
        if (sequence != nullptr)
        {
          const char *curr_start = char_finder.find(sequence,end_ptr);
          if (curr_start != nullptr)
          {
            range_start = static_cast<size_t>(curr_start - sequence);
            curr_end = char_finder.find_not(curr_start+1,end_ptr);
            if (curr_end == nullptr)
            {
              range_length = static_cast<size_t>(end_ptr - curr_start);
              curr_start = nullptr;
            } else
            {
              range_length = static_cast<size_t>(curr_end - curr_start);
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
          assert(curr_end != nullptr);
          const char *curr_start = char_finder.find(curr_end+1,end_ptr);
          if (curr_start != nullptr)
          {
            range_start = static_cast<size_t>(curr_start - sequence);
            curr_end = char_finder.find_not(curr_start+1,end_ptr);
            if (curr_end == nullptr)
            {
              range_length = static_cast<size_t>(end_ptr - curr_start);
              curr_start = nullptr;
              break;
            }
            range_length = static_cast<size_t>(curr_end - curr_start);
            break;
          } else
          {
            exhausted = true;
            break;
          }
        }
        return *this;
      }
      std::pair<size_t,size_t> operator*(void) const
      {
        return {range_start,range_length};
      }
      bool operator != (const Iterator& other) const noexcept
      {
        return exhausted != other.exhausted;
      }
      void show(void)
      {
        for (size_t idx = 0; idx <= 127; idx++)
        {
          if (char_finder.is_member(static_cast<char>(idx)))
          {
            std::cout << idx << "\t" << static_cast<char>(idx) << std::endl;
          }
        }
      }
  };
  const char *sequence;
  size_t seqlen;
  public:
    GttlCharRange(const char *_sequence,size_t seqlen) :
          sequence(_sequence),
          seqlen(seqlen) {}
  Iterator begin()
  {
    auto it = Iterator(sequence,seqlen,false);
    return it;
  }
  Iterator end()
  {
    return Iterator(nullptr,0,true);
  }
};
#endif
