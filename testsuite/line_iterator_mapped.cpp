#include <cstdlib>
#include <iostream>
#include "utilities/gttl_mmap.hpp"
#include "utilities/gttl_line_iterator.hpp"

static size_t count_lines(const char *file_part, size_t len)
{
  GttlLineIteratorFromString liter(file_part, len);
  std::string buf{};
  size_t line_num = 0;
  while (liter.next(&buf))
  {
    line_num++;
  }
  return line_num;
}

static void process_file(const char *filename)
{
  Gttlmmap<char> mapped_file(filename);
  assert(mapped_file.size() > 0);
  //constexpr const int buf_size = 1 << 14;
  const size_t parts = 2;
  size_t mid = mapped_file.size()/parts;
  assert(mid < mapped_file.size());
  const char *file_contents = mapped_file.ptr();
  const char *next_newline = static_cast<const char *>
                             (memchr(file_contents + mid,'\n',
                                     mapped_file.size() - mid));
  if (next_newline == NULL)
  {
    StrFormat msg(": second part of file beginning at offset %lu "
                  " does not contain new line character",mid);
    throw msg.str();
  }
  size_t start_part2
    = static_cast<size_t>(next_newline + 1 - file_contents);
  assert(start_part2 > 0);
  std::cout << "# part [0," << (start_part2-1)
            << "," << start_part2 << "]"
            << std::endl;
  size_t lines_part1 = count_lines(file_contents,start_part2);
  std::cout << "# lines part 1\t" << lines_part1 << std::endl;
  size_t remain = mapped_file.size() - start_part2;
  std::cout << "# part [" << start_part2 << ","
            << mapped_file.size() << "," << remain << "]" << std::endl;
  size_t lines_part2 = count_lines(file_contents + start_part2, remain);
  std::cout << "# lines part 2\t" << lines_part2 << std::endl;
  size_t lines_all = count_lines(file_contents, mapped_file.size());
  std::cout << "# lines all\t" << lines_all << std::endl;
  assert(lines_all == lines_part1 + lines_part2);
}

int main(int argc,char *argv[])
{
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " <inputfile> [inpufile ..]"
              << std::endl;
    return EXIT_FAILURE;
  }
  for (int file_idx = 1; file_idx < argc; file_idx++)
  {
    try
    {
      process_file(argv[file_idx]);
    }
    catch (std::string &msg)
    {
      std::cerr << argv[0] << ": " << msg << std::endl;
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}
