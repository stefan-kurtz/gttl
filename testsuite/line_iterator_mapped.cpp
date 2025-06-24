#include <cassert>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iostream>
#include <stdexcept>
#include "utilities/str_format.hpp"
#include "utilities/gttl_mmap.hpp"
#include "utilities/gttl_line_iterator.hpp"

static size_t count_lines(const char *file_part, size_t len)
{
  static constexpr const int buf_size = 0;
  GttlLineIterator<buf_size> liter(file_part, len);
  LineIteratorSubstring line{};
  size_t line_num = 0;
  while (liter.next(&line))
  {
    line_num++;
    line.clear();
  }
  return line_num;
}

static size_t out_fileinfo(const char *filename,const char *file_contents,
                           size_t start,size_t len)
{
  if (len == 0)
  {
    return 0;
  }
  size_t lines = count_lines(file_contents + start,len);
  std::cout << "# file\t" << filename << "\t" << start << "\t"
            << (start + len - 1) << "\t" << len << "\t" << lines << '\n';
  return lines;
}

static void process_file(const char *filename)
{
  Gttlmmap<char> mapped_file(filename);
  assert(mapped_file.size() > 0);
  const size_t parts = 2;
  size_t mid = mapped_file.size()/parts;
  assert(mid < mapped_file.size());
  const char *file_contents = mapped_file.ptr();
  const char *next_newline = static_cast<const char *>
                             (memchr(file_contents + mid,'\n',
                                     mapped_file.size() - mid));
  if (next_newline == nullptr)
  {
    StrFormat msg(": second part of file beginning at offset %zu "
                  " does not contain new line character",mid);
    throw std::runtime_error{msg.str()};
  }
  size_t start_part2
    = static_cast<size_t>(next_newline + 1 - file_contents);
  assert(start_part2 > 0);
  size_t lines1 = out_fileinfo(filename,file_contents,0,start_part2);
  size_t remain = mapped_file.size() - start_part2;
  size_t lines2 = out_fileinfo(filename,file_contents,start_part2,remain);
  size_t lines_all = count_lines(file_contents,mapped_file.size());
  std::cout << "# lines all\t" << lines_all << '\n';
  if (lines_all != lines1 + lines2)
  {
    StrFormat msg(": inconsistent number of lines: lines_all = %zu != %zu = "
                  "%zu + %zu = lines1 + lines2",lines_all,lines1+lines2,
                                                lines1,lines2);
    throw std::runtime_error{msg.str()};
  }
}

int main(int argc,char *argv[])
{
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " <inputfile> [inpufile ..]\n";
    return EXIT_FAILURE;
  }
  for (int file_idx = 1; file_idx < argc; file_idx++)
  {
    try
    {
      process_file(argv[file_idx]);
    }
    catch (const std::exception &err)
    {
      std::cerr << argv[0] << ": " << err.what() << '\n';
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}
