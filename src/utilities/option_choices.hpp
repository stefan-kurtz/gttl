#ifndef OPTION_CHOICES_HPP
#define OPTION_CHOICES_HPP
#include <algorithm>
#include <string>
#include <map>
#include <vector>
#include <iterator>
#include <sstream>

class OptionChoices
{
  private:
  const std::map<std::string, std::string> &choices_map;
  public:
  OptionChoices(const std::map<std::string, std::string> &_choices_map)
    : choices_map(_choices_map)
  { }
  [[nodiscard]] std::string help_line(void) const
  {
    std::vector<std::string> helpline_vec;
    helpline_vec.reserve(choices_map.size());
    for (auto const &[opt, helpline] : choices_map)
    {
      helpline_vec.push_back(std::string("\"").append(opt).append("\": ")
                                              .append(helpline));
    }
    const char * const delim = "; ";
    std::ostringstream help_line_os;
    std::ranges::copy(helpline_vec,
              std::ostream_iterator<std::string>(help_line_os, delim));
    return help_line_os.str();
  }
  [[nodiscard]] int choose(const std::string &choice) const
  {
    int idx = 0;
    for (auto const &[opt, helpline] : choices_map)
    {
      if (opt == choice)
      {
        return idx;
      }
      idx++;
    }
    return -1;
  }
};
#endif
