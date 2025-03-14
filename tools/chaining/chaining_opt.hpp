#ifndef UNWORDS_OPT_HPP
#define UNWORDS_OPT_HPP
#include <cstddef>
#include <string>
#include <vector>
#include <cstdbool>

class ChainingOptions
{
  private:
  std::vector<std::string> inputfiles;
  bool help_option,
       local_option,
       silent_option;

  public:
  ChainingOptions(void);
  void parse(int argc, char **argv);
  bool help_option_is_set(void) const noexcept;
  bool local_option_is_set(void) const noexcept;
  bool silent_option_is_set(void) const noexcept;
  const std::string &inputfile_get(void) const noexcept;
};
#endif
