#ifndef TAGGED_LINES_HPP
#define TAGGED_LINES_HPP
#include <ios>
#include <string>
#include <vector>
#include <fstream>

static inline std::vector<std::string> gttl_read_tagged_lines(
                                         const std::vector<std::string> &tags,
                                         const std::string &inputfile)
{
  std::ifstream infile_stream(inputfile);
  if (infile_stream.fail())
  {
    throw std::ios_base::failure(
            std::string("cannot open file ") + std::string(inputfile));
  }
  std::string line_buffer{};
  std::vector<std::string> tagged_lines{};
  while (getline(infile_stream, line_buffer))
  {
    for (auto &&tag : tags)
    {
      if (line_buffer.size() >= tag.size() &&
          line_buffer.starts_with(tag))
      {
        tagged_lines.push_back(line_buffer.substr(tag.size(),
                               line_buffer.size() - tag.size()));
        if (tagged_lines.size() == tags.size())
        {
          return tagged_lines;
        }
      }
    }
  }
  return tagged_lines;
}
#endif
