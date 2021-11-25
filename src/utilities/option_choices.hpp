#ifndef OPTION_CHOICES_HPP
#define OPTION_CHOICES_HPP
#include <string>
#include <vector>
#include <iterator>
#include <sstream>
class OptionChoices
{
  private:
  const std::vector<std::string> &choices{};
  public:
  OptionChoices(const std::vector<std::string> _choices) :
     choices(_choices)
  {
  }
  std::string help_line(void) const
  {
    const char * const delim = ", ";
    std::ostringstream help_line_os;
    std::copy(choices.begin(), choices.end(),
              std::ostream_iterator<std::string>(help_line_os, delim));
    return help_line_os.str();
  }
  int choose(const std::string &choice) const
  {
    int idx = 0;
    for (auto &&c : choices)
    {
      if (choice == c)
      {
        return idx;
      }
    }
    return -1;
  }
};
#endif
