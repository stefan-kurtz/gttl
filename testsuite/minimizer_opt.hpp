#ifndef MINIMIZER_OPT_HPP
#define MINIMIZER_OPT_HPP
#include <cstddef>
#include <string>
#include <vector>

class MinimizerOptions
{
 private:
  std::vector<std::string> inputfiles{};
  size_t qgram_length = 0,
         window_size = 1,
         number_of_threads = 1;
  int hash_bits = -1;
  bool canonical_option = false,
       at_constant_distance_option = false,
       sort_by_hash_value_option = false,
       help_option = false;
  int show_mode = 0;
  public:
  MinimizerOptions(void);
  void parse(int argc, char **argv);
  [[nodiscard]] const std::vector<std::string> &
  inputfiles_get(void) const noexcept;
  [[nodiscard]] size_t qgram_length_get(void) const noexcept;
  [[nodiscard]] size_t window_size_get(void) const noexcept;
  [[nodiscard]] size_t number_of_threads_get(void) const noexcept;
  [[nodiscard]] int hash_bits_get(void) const noexcept;
  [[nodiscard]] bool canonical_option_is_set(void) const noexcept;
  [[nodiscard]] bool at_constant_distance_option_is_set(void) const noexcept;
  [[nodiscard]] bool sort_by_hash_value_option_is_set(void) const noexcept;
  [[nodiscard]] int show_mode_get(void) const noexcept;
  [[nodiscard]] bool help_option_is_set(void) const noexcept;
};
#endif
