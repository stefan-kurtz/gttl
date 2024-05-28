#ifndef FIELDS_LINE2INDEXES_HPP
#define FIELDS_LINE2INDEXES_HPP
#include <cstddef>
#include <map>
#include <string>
#include <vector>
#include "utilities/split_string.hpp"

static inline std::vector<size_t> fields_line2indexes(
                                     const std::string &fields_line,
                                     const std::vector<std::string>
                                       &column_headers)
{
  constexpr const int skip = 2;
  std::vector<std::string> fields = gttl_split_string(fields_line,',',skip);
  size_t column = 0;
  std::map<std::string,size_t> fields_map{};
  for (auto &&f : fields)
  {
    fields_map[f] = column++;
  }
  std::vector<size_t> index_vector{};
  index_vector.reserve(column_headers.size());
  for (auto && ch : column_headers)
  {
    index_vector.push_back(fields_map[ch]);
  }
  return index_vector;
}
#endif
