#ifndef UNWORDS_OPT_HPP
#define UNWORDS_OPT_HPP
#include <cstddef>
#include <string>
#include <vector>
#include <cstdbool>

class UnwordsOptions
{
  private:
  std::vector<std::string> inputfiles;
  size_t qgram_length_max;
  bool help_option,
       ignore_reverse_complement_option,
       store_sequences_option;

  public:
  UnwordsOptions(void);
  void parse(int argc, char **argv);
  bool help_option_is_set(void) const noexcept;
  bool ignore_reverse_complement_option_is_set(void) const noexcept;
  bool store_sequences_option_is_set(void) const noexcept;
  size_t qgram_length_max_get(void) const noexcept;
  const std::vector<std::string> &inputfiles_get(void) const noexcept;
};
#endif
