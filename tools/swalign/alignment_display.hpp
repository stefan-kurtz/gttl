#ifndef ALIGNMENT_DISPLAY_HPP
#define ALIGNMENT_DISPLAY_HPP
#include <cstddef>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <cassert>
#include <vector>
#include "utilities/split_string.hpp"

class AlignmentDisplay
{
  static constexpr const int Nothing = 0,
                             Verify = 1,
                             Scoverage = (1 << 1),
                             Qcoverage = (1 << 2),
                             Identity = (1 << 3),
                             Cigar = (1 << 4),
                             Ssubstring = (1 << 5),
                             Qsubstring = (1 << 6),
                             Alignment = (1 << 7);
  int flag;
  size_t width_value;
  void set(size_t readint)
  {
    switch (readint)
    {
      case 0:
        throw std::invalid_argument("argument to option -a must be positive");
        break;
      case 1:
        flag = Verify;
        break;
      case 2:
        flag |= Scoverage;
        break;
      case 3:
        flag |= Qcoverage;
        break;
      case 4:
        flag |= Identity;
        break;
      case 5:
        flag |= Cigar;
        break;
      case 6:
        flag |= Ssubstring;
        break;
      case 7:
        flag |= Qsubstring;
        break;
      default:
        flag |= Alignment;
        width_value = readint;
        break;
    }
  }
  public:
  AlignmentDisplay(void)
    : flag(Nothing)
    , width_value(0)
  {}
  bool need_alignment(void) const
  {
    return (flag & Verify) or (flag & Identity) or
           (flag & Cigar) or (flag & Alignment);
  }
  bool need_traceback(void) const
  {
    return (flag & Identity) or (flag & Cigar) or (flag & Alignment);
  }
  bool only_verify_score(void) const
  {
    return flag == Verify;
  }
  bool identity(void) const
  {
    return flag & Identity;
  }
  bool s_coverage(void) const
  {
    return flag & Scoverage;
  }
  bool q_coverage(void) const
  {
    return flag & Qcoverage;
  }
  bool cigar(void) const
  {
    return flag & Cigar;
  }
  bool s_substring(void) const
  {
    return flag & Ssubstring;
  }
  bool q_substring(void) const
  {
    return flag & Qsubstring;
  }
  bool subject_query_alignment(void) const
  {
    return flag & Alignment;
  }
  size_t width(void) const
  {
    return width_value;
  }
  std::string to_string(void) const
  {
    if (flag & Alignment)
    {
      return std::to_string(width_value);
    }
    assert(width_value == 0);
    return std::to_string(flag);
  }
  void set_from_string(size_t min_alignment_width,const std::string &arg_opt)
  {
    const std::vector<std::string> numbers = gttl_split_string(arg_opt, '+', 1);
    size_t previous_number = 0;
    bool alignment_width_set = false;
    for (auto &s : numbers)
    {
      int readint;
      if (std::sscanf(s.c_str(),"%d",&readint) != 1 or readint < 1)
      {
        throw std::invalid_argument(
              std::string("illegal argument \"") +
              arg_opt +
              std::string("\" to option -a: + separated positive "
                          "numbers expected"));
      }
      const size_t current_number = static_cast<size_t>(readint);
      if (previous_number >= current_number)
      {
        throw std::invalid_argument(
              std::string("illegal argument ") +
              arg_opt +
              std::string(" to option -a: numbers must be strictly ordered"));
      }
      if (current_number >= min_alignment_width)
      {
        if (alignment_width_set)
        {
          throw std::invalid_argument(
                std::string("illegal argument ") +
                arg_opt +
                std::string(" to option -a: alignment width >= ") +
                std::to_string(min_alignment_width) +
                std::string(" can only be set once"));
        }
        alignment_width_set = true;
      }
      this->set(current_number);
      previous_number = current_number;
    }
  }
};
#endif
