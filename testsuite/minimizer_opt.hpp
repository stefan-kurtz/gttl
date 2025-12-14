#ifndef MINIMIZER_OPT_HPP
#define MINIMIZER_OPT_HPP
#include <cstddef>
#include <string>
#include <vector>

class MinimizerOptions
{
 private:
  std::vector<std::string> inputfiles;
  size_t qgram_length;
  size_t window_size;
  size_t number_of_threads;
  size_t max_replicates;
  int hash_bits;
  bool canonical_option;
  bool at_constant_distance_option;
  bool sort_by_hash_value_option;
  bool help_option;
  int show_mode;
  public:
  MinimizerOptions(void);
  void parse(int argc, char **argv);
  [[nodiscard]] const std::vector<std::string> &inputfiles_get(void)
                                                   const noexcept;
  [[nodiscard]] size_t qgram_length_get(void) const noexcept;
  [[nodiscard]] size_t window_size_get(void) const noexcept;
  [[nodiscard]] size_t number_of_threads_get(void) const noexcept;
  [[nodiscard]] int hash_bits_get(void) const noexcept;
  [[nodiscard]] bool canonical_option_is_set(void) const noexcept;
  [[nodiscard]] bool at_constant_distance_option_is_set(void) const noexcept;
  [[nodiscard]] bool sort_by_hash_value_option_is_set(void) const noexcept;
  [[nodiscard]] int show_mode_get(void) const noexcept;
  [[nodiscard]] size_t max_replicates_get(void) const noexcept;
  [[nodiscard]] bool help_option_is_set(void) const noexcept;
};
#endif
