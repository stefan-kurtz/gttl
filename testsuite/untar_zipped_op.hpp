#ifndef UNTAR_ZIPPED_OP_HPP
#define UNTAR_ZIPPED_OP_HPP
#include <vector>
#include <string>

class UnzippedTarOptions
{
  private:
  std::vector<std::string> inputfiles;
  bool store_option = false,
       help_option = false;
  public:
  UnzippedTarOptions(void);
  void parse(int argc, char **argv);
  const std::vector<std::string> &inputfiles_get(void) const noexcept;
  bool store_option_is_set(void) const noexcept;
  bool help_option_is_set(void) const noexcept;
};
#endif
