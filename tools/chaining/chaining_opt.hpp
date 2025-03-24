#ifndef CHAINING_OPT_HPP
#define CHAINING_OPT_HPP
#include <cstdbool>
#include <string>
#include <vector>

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
  [[nodiscard]] bool help_option_is_set(void) const noexcept;
  [[nodiscard]] bool local_option_is_set(void) const noexcept;
  [[nodiscard]] bool silent_option_is_set(void) const noexcept;
  [[nodiscard]] const std::string &inputfile_get(void) const noexcept;
};
#endif
