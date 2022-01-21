#ifndef STRING_TOKENIZER_HPP
#define STRING_TOKENIZER_HPP
#include <vector>
#include <string>
#include <cctype>

class StringTokenizer
{
  private:
    std::vector<std::string> token_list;
  public:
  StringTokenizer(const std::string &instring) :
    token_list({})
  {
    bool inword = false;
    std::string buffer{};

    for (auto it = instring.begin(); it != instring.end(); ++it)
    {
      const char cc = *it;

      if (isalpha(cc) || isdigit(cc) || cc == (int) '_')
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
