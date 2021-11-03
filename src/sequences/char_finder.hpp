#ifndef CHAR_FINDER_HPP
#define CHAR_FINDER_HPP
#include <cstdbool>

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
class MultiCharFinder
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
  template<int step,bool ref_value>
  const char *find_generic(const char *s,const char *endptr) const noexcept
  {
    for (const char *sptr = s; sptr != endptr; sptr += step)
    {
      if (in_charset[static_cast<int>(*sptr)] == ref_value)
      {
        return sptr;
      }
    }
    return nullptr;
  }
  public:
  bool is_member(char cc) const noexcept
  {
    return in_charset[static_cast<int>(cc)];
  }
  const char *find_forward(const char *s,const char *endptr) const noexcept
  {
    return find_generic<1,true>(s,endptr);
  }
  const char *find_backward(const char *s,const char *endptr) const noexcept
  {
    return find_generic<-1,true>(s,endptr);
  }
  const char *find_forward_not(const char *s,const char *endptr) const noexcept
  {
    return find_generic<1,false>(s,endptr);
  }
  const char *find_backward_not(const char *s,const char *endptr) const noexcept
  {
    return find_generic<-1,false>(s,endptr);
  }
};

template<char singlechar>
class SingleCharFinder
{
  private:
  template<int step,bool ref_value>
  const char *find_generic(const char *s,const char *endptr) const noexcept
  {
    for (const char *sptr = s; sptr != endptr; sptr += step)
    {
      if constexpr (ref_value)
      {
        if (*sptr == singlechar)
        {
          return sptr;
        }
      } else
      {
        if (*sptr != singlechar)
        {
          return sptr;
        }
      }
    }
    return nullptr;
  }
  public:
  const char *find_forward(const char *s,const char *endptr) const noexcept
  {
    return find_generic<1,true>(s,endptr);
  }
  const char *find_backward(const char *s,const char *endptr) const noexcept
  {
    return find_generic<-1,true>(s,endptr);
  }
  const char *find_forward_not(const char *s,const char *endptr) const noexcept
  {
    return find_generic<1,false>(s,endptr);
  }
  const char *find_backward_not(const char *s,const char *endptr) const noexcept
  {
    return find_generic<-1,false>(s,endptr);
  }
};
#endif
