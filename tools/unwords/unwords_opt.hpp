#ifndef UNWORDS_OPT_HPP
#define UNWORDS_OPT_HPP
#include <cstddef>
#include <string>
#include <vector>

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
  [[nodiscard]] bool help_option_is_set(void) const noexcept;
  [[nodiscard]] bool
  ignore_reverse_complement_option_is_set(void) const noexcept;
  [[nodiscard]] bool store_sequences_option_is_set(void) const noexcept;
  [[nodiscard]] size_t qgram_length_max_get(void) const noexcept;
  [[nodiscard]] const std::vector<std::string> &
  inputfiles_get(void) const noexcept;
};
#endif
