#include <cstdlib>
#include <iostream>
#include "utilities/gttl_line_generator.hpp"
#include "utilities/str_format.hpp"
#include "utilities/gttl_mmap.hpp"

static size_t count_lines(const char *file_part, size_t len)
{
  GttlLineGenerator liter(file_part, len);
  size_t line_num = 0;
  for (auto line : liter)
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
  std::cout << "# file\t" << filename << "\t"
            << start << "\t" << (start + len - 1) << "\t"
            << len << "\t" << lines << std::endl;
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
  if (next_newline == NULL)
  {
    StrFormat msg(": second part of file beginning at offset %zu "
                  " does not contain new line character",mid);
    throw msg.str();
  }
  size_t start_part2
    = static_cast<size_t>(next_newline + 1 - file_contents);
  assert(start_part2 > 0);
  size_t lines1 = out_fileinfo(filename,file_contents,0,start_part2);
  size_t remain = mapped_file.size() - start_part2;
  size_t lines2 = out_fileinfo(filename,file_contents,start_part2,remain);
  size_t lines_all = count_lines(file_contents,mapped_file.size());
  std::cout << "# lines all\t" << lines_all << std::endl;
  if (lines_all != lines1 + lines2)
  {
    StrFormat msg(": inconsistent number of lines: lines_all = %zu != %zu = "
                  "%zu + %zu = lines1 + lines2",lines_all,lines1+lines2,
                                                lines1,lines2);
    throw msg.str();
  }
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
