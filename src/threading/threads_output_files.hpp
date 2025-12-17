#ifndef THREADS_OUTPUT_FILES_HPP
#define THREADS_OUTPUT_FILES_HPP
#include <string.h> // NOLINT(modernize-deprecated-headers)
#include <stdlib.h> // NOLINT (modernize-deprecated-headers)
#include <cassert>
#include <cstddef>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <iostream>
#include <format>
#ifndef _WIN32
#include <unistd.h>
#else
#define NOMINMAX
#include <io.h>
#include "utilities/windows_mkdtemp.hpp"
#endif


class ThreadsOutputFiles
{
  bool has_threads_out_prefix;
  char *dirname_template;
  std::vector<std::string> output_filenames;
  std::vector<FILE *> output_filepointers;

  public:
  ThreadsOutputFiles(const char *program_prefix,
                     const char *threads_out_prefix,size_t num_threads)
    : has_threads_out_prefix(threads_out_prefix != nullptr)
    , dirname_template(nullptr)
    , output_filenames({})
    , output_filepointers({})
  {
    assert(num_threads > 1);
    const char *cc_threads_out_prefix;
    char separator;
    if (threads_out_prefix == nullptr)
    {
      const std::string dirname_path = std::string("./") +
                                       std::string(program_prefix) +
                                       std::string(".tmp") +
                                       std::string(".XXXXXX");
      dirname_template = strdup(dirname_path.c_str());
      cc_threads_out_prefix = mkdtemp(dirname_template);
      separator = '/';
      std::cerr << "# store output in temporary directory "
                << cc_threads_out_prefix << '\n';
    } else
    {
      cc_threads_out_prefix = threads_out_prefix;
      separator = '_';
    }
    for (size_t t_idx = 0; t_idx < num_threads; t_idx++)
    {
      const std::string fname{std::format("{}{}thread_{:02d}.tsv",
                                          cc_threads_out_prefix,
                                          separator,
                                          t_idx)};
      output_filenames.push_back(fname);
      // NOLINTNEXTLINE(misc-const-correctness)
      FILE *const out_fp = std::fopen(fname.c_str(), "w");
      if (out_fp == nullptr)
      {
        throw std::ios_base::failure(
                std::format("cannot create file \"{}\"", fname));
      }
      output_filepointers.push_back(out_fp);
    }
  }
  ~ThreadsOutputFiles(void)
  {
    for (auto *out_fp : output_filepointers)
    {
      fclose(out_fp);
    }
    if (!has_threads_out_prefix)
    {
      for (const auto &fname : output_filenames)
      {
        FILE *const in_fp     = std::fopen(fname.c_str(), "r");
        const size_t buf_size = size_t(1) << 14;
        char buf[buf_size];
        while (true)
        {
          const size_t read_bytes = fread(&buf[0], 1, buf_size, in_fp);
          std::fwrite(&buf[0],1,read_bytes,stdout);
          if (read_bytes < buf_size)
          {
            break;
          }
        }
        fclose(in_fp);
        unlink(fname.c_str());
      }
      assert(dirname_template != nullptr);
      rmdir(dirname_template);
      free(dirname_template);
    } else
    {
      for (const auto &fname : output_filenames)
      {
        printf("# output file\t%s\n",fname.c_str());
      }
    }
  }
  [[nodiscard]] FILE *filepointer(size_t t) const noexcept
  {
    assert(t < output_filepointers.size());
    return output_filepointers[t];
  }
  [[nodiscard]] const std::vector<FILE *> &
  filepointers_vector_get(void) const noexcept
  {
    return output_filepointers;
  }
  static constexpr const char *help_line =
    "specifiy prefix, say p, of files in\n"
    "       which the threads store the output. The name of the files\n"
    "       created is p_thread_tt.tsv, where tt is the 2\n"
    "       digits thread number (counting from 0) with\n"
    "       leading zeros. If this option is not specified,\n"
    "       the threads store the results in files named\n"
    "       thread_tt.tsv in a temporary directory and cat the\n"
    "       contents of the file to stdout, before removing the files and\n"
    "       the directory. It is recommended to use this option to save\n"
    "       input/output time. The names of the output\n"
    "       files are shown on stdout in lines of the form\n"
    "       # output file<tabulator>filename\n"
    "       after all files have been completly written. This\n"
    "       simplifies automatic processing of these file\n"
    "       in a pipeline without redundant specification\n"
    "       of the output files";
};
#endif
