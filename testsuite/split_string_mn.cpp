#include <cstddef>
#include <iostream>
#include <string>
#include "utilities/str_format.hpp"
#include "utilities/split_string.hpp"
#include "utilities/concatenate_strings.hpp"
#include "utilities/gttl_line_generator.hpp"

static std::pair<size_t,size_t> test_split_string(const char *inputfile,
                                                  char sep)
{
  constexpr const int buf_size = 1 << 14;
  GttlLineGenerator<buf_size> gttl_li(inputfile);
  size_t line_count = 0, column_count = 0;
  for (auto &buffer : gttl_li)
  {
    if (buffer.size() > 0 and buffer[0] != '#')
    {
      std::vector<std::string> vec = gttl_split_string(buffer, sep);
      column_count += vec.size();
      const std::string sep_string{sep};
      std::string line_from_vec = gttl_concatenate_strings(vec.begin(),
                                                           vec.end(),
                                                           sep_string);
      if (buffer != line_from_vec)
      {
        StrFormat msg(": '%s' != '%s'",buffer.c_str(),line_from_vec.c_str());
        throw msg.str();
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
    std::cerr << "Usage: " << argv[0] << " <sep> <filename>" << std::endl;
    return EXIT_FAILURE;
  }
  const char *inputfile = argv[2];
  if (std::string(argv[1]).size() > 1)
  {
    std::cerr << argv[0] << ": seperator must be single character" << std::endl;
    return EXIT_FAILURE;
  }
  try
  {
    size_t line_count, column_count;
    std::tie(line_count,column_count) = test_split_string(inputfile,argv[1][0]);
    std::cout << "processed " << column_count << " columns in " << line_count
              << " lines" << std::endl;
  }
  catch (const std::string &msg)
  {
    std::cerr << argv[0] << msg << " " << inputfile << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
