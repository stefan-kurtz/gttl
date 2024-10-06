#ifndef NTCARD_OPT_HPP
#define NTCARD_OPT_HPP
#include <cstddef>
#include <string>

class NtcardOptions
{
  private:
    std::string inputfile;
    size_t qgram_length;
    size_t s;
    size_t r;
    bool show_f_option;
    bool fast_option;
    bool binary_option;
    bool handle_wildcard_like_A;
    bool help_option;
    size_t num_threads;

  public:
    NtcardOptions(void);
    void parse(int argc, char **argv);
    const std::string &inputfile_get(void) const noexcept;
    size_t qgram_length_get(void) const noexcept;
    size_t s_get(void) const noexcept;
    size_t r_get(void) const noexcept;
    bool show_f_option_is_set(void) const noexcept;
    bool fast_option_is_set(void) const noexcept;
    bool binary_option_is_set(void) const noexcept;
    bool handle_wildcard_like_A_is_set(void) const noexcept;
    bool help_option_is_set(void) const noexcept;
    size_t num_threads_get(void) const noexcept;
};
#endif
