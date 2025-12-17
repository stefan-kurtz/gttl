#ifndef STRING_TOKENIZER_HPP
#define STRING_TOKENIZER_HPP
#include <vector>
#include <string>
#include <cctype>

class StringTokenizer
{
  std::vector<std::string> token_list;
  public:
  StringTokenizer(const std::string &instring)
  {
    bool inword = false;
    std::string buffer;

    for (const char cc : instring)
    {
      /* split string in maximum length substrings of alphabetic/digit/_/.
         and store each substring in vector */
      static constexpr const int underscore = static_cast<int>('_');
      static constexpr const int dot = static_cast<int>('.');

      if (isalpha(cc) || isdigit(cc) || cc == underscore || cc == dot)
      {
        buffer.push_back((char) cc);
        inword = true;
      } else
      {
        if (inword)
        {
          if (not buffer.empty())
          {
            token_list.emplace_back(buffer);
            buffer.clear();
          }
          inword = false;
        }
      }
    }
    if (inword and not buffer.empty())
    {
      token_list.emplace_back(buffer);
    }
  }
  [[nodiscard]] const std::vector<std::string> &
  token_list_get(void) const noexcept
  {
    return token_list;
  }
};
#endif
