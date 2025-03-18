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

    for (auto it = instring.begin(); it != instring.end(); ++it)
    {
      /* split string in maximum length substrings of alphabetic/digit/_/.
         and store each substring in vector */
      const char cc = *it;
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
          if (buffer.size() > 0)
          {
            token_list.push_back(std::string(buffer));
            buffer.clear();
          }
          inword = false;
        }
      }
    }
    if (inword && buffer.size() > 0)
    {
      token_list.push_back(std::string(buffer));
    }
  }
  const std::vector<std::string> &token_list_get(void) const noexcept
  {
    return token_list;
  }
};
#endif
