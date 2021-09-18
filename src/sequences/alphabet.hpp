#ifndef ALPHABET_HPP
#define ALPHABET_HPP
#include <climits>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <array>

template <uint8_t undef,const char *str>
static constexpr uint8_t find_index(char cc,std::size_t i,std::size_t j)
{
  return str[i] == '\0'
           ? undef
           : (str[i] == cc ? static_cast<uint8_t>(j)
                           : (str[i] == '|' ? find_index<undef,str>(cc,i+1,j+1)
                                            : find_index<undef,str>(cc,i+1,j)));
}

template <uint8_t undef,const char *str>
static constexpr uint8_t find_index (char cc)
{
  return find_index<undef,str>(cc,0,0);
}

template <const char *str>
static constexpr size_t find_size(std::size_t i,std::size_t j)
{
  return str[i] == '\0' ? j : (str[i] == '|' ? find_size<str>(i+1,j+1)
                                             : find_size<str>(i+1,j));
}

template <const char *str>
static constexpr size_t find_size(void)
{
  return 1 + find_size<str>(0,0);
}

template <size_t _size,const char *str>
static constexpr void fill_characters (std::array<char,_size> &_chars,
                                       std::size_t i,std::size_t j,
                                       std::size_t k)
{
  if (str[i] != '\0')
  {
    if (j == 0)
    {
      _chars[k] = str[i];
      fill_characters<_size,str> (_chars,i+1,j+1,k+1);
    } else
    {
      if (str[i] == '|')
      {
        fill_characters<_size,str> (_chars,i+1,0,k);
      } else
      {
        fill_characters<_size,str> (_chars,i+1,j+1,k);
      }
    }
  }
}

template <size_t _size,const char *str>
static constexpr std::array<char,_size> fill_characters (void)
{
  std::array<char,_size> _chars{};
  fill_characters<_size,str> (_chars,0,0,0);
  return _chars;
}

template <const char *char_spec,uint8_t _undefined_rank>
class GttlAlphabet
{
  static constexpr const uint8_t _symbolmap[UCHAR_MAX+1] =
  {
    /* 0 */ _undefined_rank,
    /* 1 */ _undefined_rank,
    /* 2 */ _undefined_rank,
    /* 3 */ _undefined_rank,
    /* 4 */ _undefined_rank,
    /* 5 */ _undefined_rank,
    /* 6 */ _undefined_rank,
    /* 7 */ _undefined_rank,
    /* 8 */ _undefined_rank,
    /* 9 */ _undefined_rank,
    /* 10 */ _undefined_rank,
    /* 11 */ _undefined_rank,
    /* 12 */ _undefined_rank,
    /* 13 */ _undefined_rank,
    /* 14 */ _undefined_rank,
    /* 15 */ _undefined_rank,
    /* 16 */ _undefined_rank,
    /* 17 */ _undefined_rank,
    /* 18 */ _undefined_rank,
    /* 19 */ _undefined_rank,
    /* 20 */ _undefined_rank,
    /* 21 */ _undefined_rank,
    /* 22 */ _undefined_rank,
    /* 23 */ _undefined_rank,
    /* 24 */ _undefined_rank,
    /* 25 */ _undefined_rank,
    /* 26 */ _undefined_rank,
    /* 27 */ _undefined_rank,
    /* 28 */ _undefined_rank,
    /* 29 */ _undefined_rank,
    /* 30 */ _undefined_rank,
    /* 31 */ _undefined_rank,
    /* 32 */ _undefined_rank,
    /* 33 */ find_index<_undefined_rank,char_spec>(33),
    /* 34 */ find_index<_undefined_rank,char_spec>(34),
    /* 35 */ find_index<_undefined_rank,char_spec>(35),
    /* 36 */ find_index<_undefined_rank,char_spec>(36),
    /* 37 */ find_index<_undefined_rank,char_spec>(37),
    /* 38 */ find_index<_undefined_rank,char_spec>(38),
    /* 39 */ find_index<_undefined_rank,char_spec>(39),
    /* 40 */ find_index<_undefined_rank,char_spec>(40),
    /* 41 */ find_index<_undefined_rank,char_spec>(41),
    /* 42 */ find_index<_undefined_rank,char_spec>(42),
    /* 43 */ find_index<_undefined_rank,char_spec>(43),
    /* 44 */ find_index<_undefined_rank,char_spec>(44),
    /* 45 */ find_index<_undefined_rank,char_spec>(45),
    /* 46 */ find_index<_undefined_rank,char_spec>(46),
    /* 47 */ find_index<_undefined_rank,char_spec>(47),
    /* 48 */ find_index<_undefined_rank,char_spec>(48),
    /* 49 */ find_index<_undefined_rank,char_spec>(49),
    /* 50 */ find_index<_undefined_rank,char_spec>(50),
    /* 51 */ find_index<_undefined_rank,char_spec>(51),
    /* 52 */ find_index<_undefined_rank,char_spec>(52),
    /* 53 */ find_index<_undefined_rank,char_spec>(53),
    /* 54 */ find_index<_undefined_rank,char_spec>(54),
    /* 55 */ find_index<_undefined_rank,char_spec>(55),
    /* 56 */ find_index<_undefined_rank,char_spec>(56),
    /* 57 */ find_index<_undefined_rank,char_spec>(57),
    /* 58 */ find_index<_undefined_rank,char_spec>(58),
    /* 59 */ find_index<_undefined_rank,char_spec>(59),
    /* 60 */ find_index<_undefined_rank,char_spec>(60),
    /* 61 */ find_index<_undefined_rank,char_spec>(61),
    /* 62 */ find_index<_undefined_rank,char_spec>(62),
    /* 63 */ find_index<_undefined_rank,char_spec>(63),
    /* 64 */ find_index<_undefined_rank,char_spec>(64),
    /* 65 */ find_index<_undefined_rank,char_spec>(65),
    /* 66 */ find_index<_undefined_rank,char_spec>(66),
    /* 67 */ find_index<_undefined_rank,char_spec>(67),
    /* 68 */ find_index<_undefined_rank,char_spec>(68),
    /* 69 */ find_index<_undefined_rank,char_spec>(69),
    /* 70 */ find_index<_undefined_rank,char_spec>(70),
    /* 71 */ find_index<_undefined_rank,char_spec>(71),
    /* 72 */ find_index<_undefined_rank,char_spec>(72),
    /* 73 */ find_index<_undefined_rank,char_spec>(73),
    /* 74 */ find_index<_undefined_rank,char_spec>(74),
    /* 75 */ find_index<_undefined_rank,char_spec>(75),
    /* 76 */ find_index<_undefined_rank,char_spec>(76),
    /* 77 */ find_index<_undefined_rank,char_spec>(77),
    /* 78 */ find_index<_undefined_rank,char_spec>(78),
    /* 79 */ find_index<_undefined_rank,char_spec>(79),
    /* 80 */ find_index<_undefined_rank,char_spec>(80),
    /* 81 */ find_index<_undefined_rank,char_spec>(81),
    /* 82 */ find_index<_undefined_rank,char_spec>(82),
    /* 83 */ find_index<_undefined_rank,char_spec>(83),
    /* 84 */ find_index<_undefined_rank,char_spec>(84),
    /* 85 */ find_index<_undefined_rank,char_spec>(85),
    /* 86 */ find_index<_undefined_rank,char_spec>(86),
    /* 87 */ find_index<_undefined_rank,char_spec>(87),
    /* 88 */ find_index<_undefined_rank,char_spec>(88),
    /* 89 */ find_index<_undefined_rank,char_spec>(89),
    /* 90 */ find_index<_undefined_rank,char_spec>(90),
    /* 91 */ find_index<_undefined_rank,char_spec>(91),
    /* 92 */ find_index<_undefined_rank,char_spec>(92),
    /* 93 */ find_index<_undefined_rank,char_spec>(93),
    /* 94 */ find_index<_undefined_rank,char_spec>(94),
    /* 95 */ find_index<_undefined_rank,char_spec>(95),
    /* 96 */ find_index<_undefined_rank,char_spec>(96),
    /* 97 */ find_index<_undefined_rank,char_spec>(97),
    /* 98 */ find_index<_undefined_rank,char_spec>(98),
    /* 99 */ find_index<_undefined_rank,char_spec>(99),
    /* 100 */ find_index<_undefined_rank,char_spec>(100),
    /* 101 */ find_index<_undefined_rank,char_spec>(101),
    /* 102 */ find_index<_undefined_rank,char_spec>(102),
    /* 103 */ find_index<_undefined_rank,char_spec>(103),
    /* 104 */ find_index<_undefined_rank,char_spec>(104),
    /* 105 */ find_index<_undefined_rank,char_spec>(105),
    /* 106 */ find_index<_undefined_rank,char_spec>(106),
    /* 107 */ find_index<_undefined_rank,char_spec>(107),
    /* 108 */ find_index<_undefined_rank,char_spec>(108),
    /* 109 */ find_index<_undefined_rank,char_spec>(109),
    /* 110 */ find_index<_undefined_rank,char_spec>(110),
    /* 111 */ find_index<_undefined_rank,char_spec>(111),
    /* 112 */ find_index<_undefined_rank,char_spec>(112),
    /* 113 */ find_index<_undefined_rank,char_spec>(113),
    /* 114 */ find_index<_undefined_rank,char_spec>(114),
    /* 115 */ find_index<_undefined_rank,char_spec>(115),
    /* 116 */ find_index<_undefined_rank,char_spec>(116),
    /* 117 */ find_index<_undefined_rank,char_spec>(117),
    /* 118 */ find_index<_undefined_rank,char_spec>(118),
    /* 119 */ find_index<_undefined_rank,char_spec>(119),
    /* 120 */ find_index<_undefined_rank,char_spec>(120),
    /* 121 */ find_index<_undefined_rank,char_spec>(121),
    /* 122 */ find_index<_undefined_rank,char_spec>(122),
    /* 123 */ find_index<_undefined_rank,char_spec>(123),
    /* 124 */ _undefined_rank,
    /* 125 */ find_index<_undefined_rank,char_spec>(125),
    /* 126 */ find_index<_undefined_rank,char_spec>(126),
    /* 127 */ _undefined_rank,
    /* 128 */ _undefined_rank,
    /* 129 */ _undefined_rank,
    /* 130 */ _undefined_rank,
    /* 131 */ _undefined_rank,
    /* 132 */ _undefined_rank,
    /* 133 */ _undefined_rank,
    /* 134 */ _undefined_rank,
    /* 135 */ _undefined_rank,
    /* 136 */ _undefined_rank,
    /* 137 */ _undefined_rank,
    /* 138 */ _undefined_rank,
    /* 139 */ _undefined_rank,
    /* 140 */ _undefined_rank,
    /* 141 */ _undefined_rank,
    /* 142 */ _undefined_rank,
    /* 143 */ _undefined_rank,
    /* 144 */ _undefined_rank,
    /* 145 */ _undefined_rank,
    /* 146 */ _undefined_rank,
    /* 147 */ _undefined_rank,
    /* 148 */ _undefined_rank,
    /* 149 */ _undefined_rank,
    /* 150 */ _undefined_rank,
    /* 151 */ _undefined_rank,
    /* 152 */ _undefined_rank,
    /* 153 */ _undefined_rank,
    /* 154 */ _undefined_rank,
    /* 155 */ _undefined_rank,
    /* 156 */ _undefined_rank,
    /* 157 */ _undefined_rank,
    /* 158 */ _undefined_rank,
    /* 159 */ _undefined_rank,
    /* 160 */ _undefined_rank,
    /* 161 */ _undefined_rank,
    /* 162 */ _undefined_rank,
    /* 163 */ _undefined_rank,
    /* 164 */ _undefined_rank,
    /* 165 */ _undefined_rank,
    /* 166 */ _undefined_rank,
    /* 167 */ _undefined_rank,
    /* 168 */ _undefined_rank,
    /* 169 */ _undefined_rank,
    /* 170 */ _undefined_rank,
    /* 171 */ _undefined_rank,
    /* 172 */ _undefined_rank,
    /* 173 */ _undefined_rank,
    /* 174 */ _undefined_rank,
    /* 175 */ _undefined_rank,
    /* 176 */ _undefined_rank,
    /* 177 */ _undefined_rank,
    /* 178 */ _undefined_rank,
    /* 179 */ _undefined_rank,
    /* 180 */ _undefined_rank,
    /* 181 */ _undefined_rank,
    /* 182 */ _undefined_rank,
    /* 183 */ _undefined_rank,
    /* 184 */ _undefined_rank,
    /* 185 */ _undefined_rank,
    /* 186 */ _undefined_rank,
    /* 187 */ _undefined_rank,
    /* 188 */ _undefined_rank,
    /* 189 */ _undefined_rank,
    /* 190 */ _undefined_rank,
    /* 191 */ _undefined_rank,
    /* 192 */ _undefined_rank,
    /* 193 */ _undefined_rank,
    /* 194 */ _undefined_rank,
    /* 195 */ _undefined_rank,
    /* 196 */ _undefined_rank,
    /* 197 */ _undefined_rank,
    /* 198 */ _undefined_rank,
    /* 199 */ _undefined_rank,
    /* 200 */ _undefined_rank,
    /* 201 */ _undefined_rank,
    /* 202 */ _undefined_rank,
    /* 203 */ _undefined_rank,
    /* 204 */ _undefined_rank,
    /* 205 */ _undefined_rank,
    /* 206 */ _undefined_rank,
    /* 207 */ _undefined_rank,
    /* 208 */ _undefined_rank,
    /* 209 */ _undefined_rank,
    /* 210 */ _undefined_rank,
    /* 211 */ _undefined_rank,
    /* 212 */ _undefined_rank,
    /* 213 */ _undefined_rank,
    /* 214 */ _undefined_rank,
    /* 215 */ _undefined_rank,
    /* 216 */ _undefined_rank,
    /* 217 */ _undefined_rank,
    /* 218 */ _undefined_rank,
    /* 219 */ _undefined_rank,
    /* 220 */ _undefined_rank,
    /* 221 */ _undefined_rank,
    /* 222 */ _undefined_rank,
    /* 223 */ _undefined_rank,
    /* 224 */ _undefined_rank,
    /* 225 */ _undefined_rank,
    /* 226 */ _undefined_rank,
    /* 227 */ _undefined_rank,
    /* 228 */ _undefined_rank,
    /* 229 */ _undefined_rank,
    /* 230 */ _undefined_rank,
    /* 231 */ _undefined_rank,
    /* 232 */ _undefined_rank,
    /* 233 */ _undefined_rank,
    /* 234 */ _undefined_rank,
    /* 235 */ _undefined_rank,
    /* 236 */ _undefined_rank,
    /* 237 */ _undefined_rank,
    /* 238 */ _undefined_rank,
    /* 239 */ _undefined_rank,
    /* 240 */ _undefined_rank,
    /* 241 */ _undefined_rank,
    /* 242 */ _undefined_rank,
    /* 243 */ _undefined_rank,
    /* 244 */ _undefined_rank,
    /* 245 */ _undefined_rank,
    /* 246 */ _undefined_rank,
    /* 247 */ _undefined_rank,
    /* 248 */ _undefined_rank,
    /* 249 */ _undefined_rank,
    /* 250 */ _undefined_rank,
    /* 251 */ _undefined_rank,
    /* 252 */ _undefined_rank,
    /* 253 */ _undefined_rank,
    /* 254 */ _undefined_rank,
    /* 255 */ _undefined_rank
  };
  static constexpr size_t _size = find_size<char_spec>();
  static constexpr std::array<char,_size> characters
    = fill_characters<_size,char_spec>();
  public:
    constexpr GttlAlphabet(void) { }
    constexpr uint8_t undefined_rank(void) const noexcept
    {
      return _undefined_rank;
    }
    constexpr size_t size(void) const noexcept { return _size; }
    constexpr uint8_t char_to_rank(char cc) const noexcept {
      return _symbolmap[static_cast<int>(cc)];
    }
    constexpr char rank_to_char(uint8_t r) const noexcept {
      return characters[r];
    }
    void pretty_print(void)
    {
      std::cout << "# alphabet size\t" << _size << std::endl;
      std::cout << "# undefined_rank\t" << static_cast<int>(_undefined_rank)
                << std::endl;
      for (const char *s = char_spec; *s != '\0'; s++)
      {
        if (*s != '|')
        {
          uint8_t r = this->char_to_rank(*s);
          std::cout << *s << "\t" << static_cast<int>(r) << std::endl;
        }
      }
      for (size_t idx = 0; idx < _size; idx++)
      {
        std::cout << idx << "\t" << rank_to_char(idx) << std::endl;
      }
    }
    constexpr const char *characters_get(void) const noexcept
    {
      return characters.data();
    }

    void char_to_rank_in_place(char *char_seq,size_t len) const noexcept
    {
      uint8_t *byte_seq = reinterpret_cast<uint8_t *>(char_seq);
      for (size_t idx = 0; idx < len; idx++)
      {
        byte_seq[idx] = char_to_rank(char_seq[idx]);
      }
    }
};

namespace alphabet
{
  static constexpr const char nucleotides_upper_lower[] = "Aa|Cc|Gg|TtUu";
  using GttlAlphabet_UL_4 = GttlAlphabet<nucleotides_upper_lower,4>;
  using GttlAlphabet_UL_0 = GttlAlphabet<nucleotides_upper_lower,0>;
  static constexpr const char nucleotides_upper[] = "A|C|G|TU";
  using GttlAlphabet_U_4 = GttlAlphabet<nucleotides_upper_lower,4>;
  static constexpr const char amino_acids[]
    = "A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y";
  using GttlAlphabet_20 = GttlAlphabet<amino_acids,20>;
}
#endif
