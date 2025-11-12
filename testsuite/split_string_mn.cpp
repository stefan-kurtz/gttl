#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include "utilities/str_format.hpp"
#include "utilities/split_string.hpp"
#include "utilities/concatenate_strings.hpp"
#include "utilities/gttl_line_generator.hpp"

static std::pair<size_t,size_t> test_split_string(const char *inputfile,
                                                  char sep)
{
  constexpr const int buf_size = 1 << 14;
  GttlLineGenerator<buf_size> gttl_lg(inputfile);
  size_t line_count = 0;
  size_t column_count = 0;
  for (const auto &buffer : gttl_lg)
  {
    if (buffer.size() > 0 and buffer[0] != '#')
    {
      std::vector<std::string> vec = gttl_split_string(buffer, sep);
      column_count += vec.size();
      const std::string sep_string{sep};
      const std::string line_from_vec = gttl_concatenate_strings(
                                   vec.begin(), vec.end(), sep_string);
      if (buffer != line_from_vec)
      {
        const StrFormat msg(": '%s' != '%s'",
                            buffer.c_str(),
                            line_from_vec.c_str());
        throw std::runtime_error{msg.str()};
      }
    }
    line_count++;
  }
  return std::make_pair(line_count,column_count);
}

int main(int argc,char *argv[])
{
  if (argc != 3)
  {
    std::cerr << "Usage: " << argv[0] << " <sep> <filename>\n";
    return EXIT_FAILURE;
  }
  const char *const inputfile = argv[2];
  if (std::string(argv[1]).size() > 1)
  {
    std::cerr << argv[0] << ": seperator must be single character\n";
    return EXIT_FAILURE;
  }
  try
  {
    size_t line_count = 0;
    size_t column_count = 0;
    std::tie(line_count,column_count) = test_split_string(inputfile,argv[1][0]);
    std::cout << "processed " << column_count << " columns in " << line_count
              << " lines\n";
  }
  catch (const std::exception &err)
  {
    std::cerr << argv[0] << err.what() << " " << inputfile << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
