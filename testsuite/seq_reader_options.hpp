#ifndef SEQ_READER_OPTIONS_HPP
#define SEQ_READER_OPTIONS_HPP
#include <cstddef>
#include <cstdbool>
#include <string>
#include <vector>

class SeqReaderOptions
{
 private:
  std::vector<std::string> inputfiles{};
  bool help_option = false,
       statistics_option = false,
       echo_option = false,
       fasta_output_option = false,
       mapped_option = false;
  size_t num_threads = 0, split_size = 0;

 public:
  SeqReaderOptions(void);

  void parse(int argc, char **argv);
  bool help_option_is_set(void) const noexcept;
  bool statistics_option_is_set(void) const noexcept;
  bool echo_option_is_set(void) const noexcept;
  bool fasta_output_option_is_set(void) const noexcept;
  bool mapped_option_is_set(void) const noexcept;
  size_t split_size_get(void) const noexcept;
  size_t num_threads_get(void) const noexcept;
  const std::vector<std::string> &inputfiles_get(void) const noexcept;
};
#endif
