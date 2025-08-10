#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include "utilities/find_lit_string.hpp"
#include "utilities/multiple_options.hpp"

static constexpr const GttlLitStringInitializerList display_option_args
{
  "k_mer",
  "k_mer_counts",
  "log2FC",
  "p_value",
  "p_value_histogram"
};

using DisplayOptions = MultipleOptions<display_option_args>;

int main(int argc, char *argv[])
{
  DisplayOptions display_options("");
  if (argc != 2)
  {
    std::cerr << argv[0] << (": missing argument specifying + separated values "
                             "in single string, possible values are: ")
                         << display_options.help_string()
                         << '\n';
    return EXIT_FAILURE;
  }
  const std::string argstring(argv[1]);
  try
  {
    display_options.set_flags(argstring);
    std::cout << argv[0] << ": the following values are set: "
                         <<  display_options.to_string() << '\n';
    if (display_options.is_set<"p_value">())
    {
      std::cout << "p_value is set" << '\n';
    }
  }
  catch (const std::invalid_argument &err)
  {
    std::cerr << argv[0] << err.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
