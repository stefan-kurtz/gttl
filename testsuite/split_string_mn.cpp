#include <cstddef>
#include <iostream>
#include <string>
#include "utilities/split_string.hpp"
#include "utilities/gttl_line_iterator.hpp"

static std::string concatenate(const std::vector<std::string> &vec,char sep)
{
  std::string string_from_vec{};
  for (auto &s : vec)
  {
    if (string_from_vec.size() > 0)
    {
      string_from_vec += sep;
    }
    string_from_vec += s;
  }
  string_from_vec += '\n';
  return string_from_vec;
}


static std::pair<size_t,size_t> process_inputfile(const char *inputfile,
                                                  char sep)
{
  constexpr const int buf_size = 1 << 14;
  GttlLineIterator<buf_size> gttl_li(inputfile);
  std::string buffer{};
  size_t line_count = 0, column_count = 0;
  while (gttl_li.next(&buffer))
  {
    if (buffer.size() > 0 and buffer[0] != '#')
    {
      std::vector<std::string> vec = gttl_split_string(buffer, sep);
      column_count += vec.size();
      std::string line_from_vec = concatenate(vec,sep);
      if (buffer != line_from_vec)
      {
        StrFormat msg("'%s' != '%s'",buffer.c_str(),line_from_vec.c_str());
        throw msg.str();
      }
    }
    buffer.clear();
    line_count++;
  }
  return {line_count,column_count};
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
    std::tie(line_count,column_count) = process_inputfile(inputfile,argv[1][0]);
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
