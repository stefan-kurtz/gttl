#include <iostream>
#include <stdexcept>
#include <vector>
#include "utilities/constexpr_for.hpp"
#include "utilities/gttl_line_generator.hpp"
#include "utilities/split_string.hpp"
#include "stored_match.hpp"
#include "chaining.hpp"
#include "chaining_opt.hpp"

int main(int argc, char *argv[])
{
  ChainingOptions options;
  try
  {
    options.parse(argc, argv);
  }
  catch (const std::invalid_argument &e)
  {
    std::cerr << argv[0] << e.what() << '\n';
    return EXIT_FAILURE;
  }
  if (options.help_option_is_set())
  {
    return EXIT_SUCCESS;
  }

  constexpr const int buf_size = 1U << 14U;
  const std::string inputfile = options.inputfile_get();
  const bool local_option = options.local_option_is_set();
  const bool silent_option = options.silent_option_is_set();
  std::vector<GttlStoredMatch> matches{};
  try
  {
    GttlLineGenerator<buf_size> gttl_lg(inputfile.c_str());
    for(const auto& line : gttl_lg)
    {
      std::vector<std::string> vec = gttl_split_string(line, ' ');
      if (vec.size() != 5)
      {
        throw std::runtime_error(": line has " + std::to_string(vec.size())
                                 + " columns, but 5 are expected");
      }
      matches.emplace_back(std::stoul(vec[0]),
                           std::stoul(vec[1])-std::stoul(vec[0])+1,
                           std::stoul(vec[2]),
                           std::stoul(vec[3])-std::stoul(vec[2])+1,
                           0,
                           std::stoul(vec[4]));
    }
  }
  catch (std::string &msg)
  {
    std::cerr << argv[0] << ": " << msg << '\n';
    return EXIT_FAILURE;
  }

  constexpr_for<0,1+1,1>([&](auto compile_time_local_option)
  {
    if (compile_time_local_option == local_option)
    {
      Chain<GttlStoredMatch, compile_time_local_option> chaining(matches);
      std::cout << "# chain: length " << chaining.size()
                << " score " << chaining.score() << '\n';
      if (not silent_option)
      {
        for (auto idx : chaining)
        {
          std::cout << idx << " " << matches[idx].to_string() << '\n';
        }
      }
    }
  });
  return EXIT_SUCCESS;
}
