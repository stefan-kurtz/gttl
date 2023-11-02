#ifndef MINIMIZER_OPT_HPP
#define MINIMIZER_OPT_HPP
#include <cstddef>
#include <string>
#include <vector>

class MinimizerOptions
{
 private:
  std::vector<std::string> inputfiles;
  size_t qgram_length,
         window_size,
         number_of_threads;
  int hash_bits;
  bool canonical_option,
       at_constant_distance_option,
       sort_by_hash_value_option,
       help_option;
  public:
  MinimizerOptions(void);
  void parse(int argc, char **argv);
  const std::vector<std::string> &inputfiles_get(void) const noexcept;
  size_t qgram_length_get(void) const noexcept;
  size_t window_size_get(void) const noexcept;
  size_t number_of_threads_get(void) const noexcept;
  int hash_bits_get(void) const noexcept;
  bool canonical_option_is_set(void) const noexcept;
  bool at_constant_distance_option_is_set(void) const noexcept;
  bool sort_by_hash_value_option_is_set(void) const noexcept;
  bool help_option_is_set(void) const noexcept;
};
#endif
