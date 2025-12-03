#ifndef DUPLICATED_FILENAMES_HPP
#define DUPLICATED_FILENAMES_HPP

#include <stdlib.h> // NOLINT(modernize-deprecated-headers)
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>
#include <format>

static inline bool gttl_files_with_same_path(const std::string &path0,
                                             const std::string &path1)
{
  char *const fullpath_0 = realpath(path0.c_str(), nullptr);
  char *const fullpath_1 = realpath(path1.c_str(), nullptr);
  const bool same_files = std::strcmp(fullpath_0,
                                      fullpath_1) == 0;
  free(fullpath_0);
  free(fullpath_1);
  return same_files;
}

static inline void gttl_duplicated_filenames(const std::vector<std::string>
                                                &input_files)
{
  for (size_t i = 0; i < input_files.size(); i++)
  {
    for (size_t j = i+1; j < input_files.size(); j++)
    {
      if (gttl_files_with_same_path(input_files[i],input_files[j]))
      {
        throw std::runtime_error(
                std::format("{} appears twice in list of input files; if you "
                            "want to compare each sequence of a file against "
                            "each other (e.g. to search for repeats), then run "
                            "the program only with this single filename",
                            input_files[i]));
      }
    }
  }
}
#endif
