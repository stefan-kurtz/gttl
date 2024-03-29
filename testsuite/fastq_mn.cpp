/*
  Copyright (c) 2021 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2021 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/
#include <cstddef>
#include <iostream>
#include <fstream>
#include <cstdbool>
#include <typeinfo>
#include <cassert>
#include <algorithm>
#include <thread>
#include <cstdio>
#include "utilities/str_format.hpp"
#include "utilities/basename.hpp"
#include "utilities/gttl_mmap.hpp"
#include "sequences/gttl_fastq_iterator.hpp"
#include "sequences/fastq_parts.hpp"
#include "seq_reader_options.hpp"

static void fastq_split_writer(size_t split_size,
                               const std::string &inputfilename)
{
  constexpr const int buf_size = 1 << 14;
  GttlLineIterator<buf_size> line_iterator(inputfilename.c_str());
  GttlFastQIterator<GttlLineIterator<buf_size>> fastq_it(line_iterator);
  int file_number = 0;
  auto it = fastq_it.begin();
  bool exhausted = false;
  GttlBasename inputfile_base(inputfilename.c_str());
  while (!exhausted)
  {
    if (it == fastq_it.end())
    {
      break;
    }
    StrFormat outfilename("%s_%02d",inputfile_base.str(),file_number++);
    std::ofstream out_stream;
    out_stream.open(outfilename.str());
    for (size_t idx = 0; idx < split_size; idx++)
    {
      if (it == fastq_it.end())
      {
        exhausted = true;
        break;
      }
      const std::string_view &sequence = (*it).sequence_get();
      const std::string_view &header = (*it).header_get();
      const std::string_view &quality = (*it).quality_get();
      out_stream << header << std::endl
                 << sequence << std::endl
                 << "+" << std::endl
                 << quality << std::endl;
      ++it;
    }
    out_stream.close();
  }
}

template<class FastQIterator>
static void process_fastq_iter(bool statistics,
                               bool echo,
                               bool fasta_output,
                               FastQIterator &fastq_it)
{
  size_t seqnum = 0, total_length = 0;
  for (auto &&fastq_entry : fastq_it)
  {
    const std::string_view &sequence = fastq_entry.sequence_get();
    if (echo)
    {
      const std::string_view &header = fastq_entry.header_get();
      const std::string_view &quality = fastq_entry.quality_get();
      std::cout << header << std::endl
                << sequence << std::endl
                << "+" << std::endl
                << quality << std::endl;
    } else
    {
      if (fasta_output)
      {
        const std::string_view &header = fastq_entry.header_get();
        std::cout << ">" << header.substr(1) << std::endl
                  << sequence << std::endl;
      }
    }
    total_length += sequence.size();
    seqnum++;
  }
  if (statistics)
  {
    std::cout << "# number of sequences\t" << seqnum << std::endl;
    std::cout << "# total length\t" << total_length << std::endl;
    std::cout << "# mean length\t" << total_length/seqnum << std::endl;
  }
}

static void process_single_file_streamed(bool statistics,
                                         bool echo,
                                         bool fasta_output,
                                         const std::string &inputfilename)
{
  constexpr const int buf_size = 1 << 14;
  GttlLineIterator<buf_size> line_iterator(inputfilename.c_str());
  using FastQIterator = GttlFastQIterator<GttlLineIterator<buf_size>>;
  FastQIterator fastq_it(line_iterator);
  process_fastq_iter<FastQIterator>(statistics,echo,fasta_output,fastq_it);
}

static void process_single_file_mapped(bool statistics,
                                       bool echo,
                                       bool fasta_output,
                                       const std::string &inputfilename)
{
  constexpr const int buf_size = 0;
  Gttlmmap<char> mapped_file(inputfilename.c_str());
  GttlLineIterator<buf_size> line_iterator(mapped_file.ptr(),
                                           mapped_file.size());
  using FastQIterator = GttlFastQIterator<GttlLineIterator<buf_size>>;
  FastQIterator fastq_it(line_iterator);
  process_fastq_iter<FastQIterator>(statistics,echo,fasta_output,fastq_it);
}

static void process_paired_files(bool statistics,
                                 bool fasta_output,
                                 const std::string &filename0,
                                 const std::string &filename1)
{
  constexpr const int buf_size = 1 << 14;
  GttlLineIterator<buf_size> line_iterator0(filename0.c_str()),
                             line_iterator1(filename1.c_str());
  using FastQIterator = GttlFastQIterator<GttlLineIterator<buf_size>>;
  FastQIterator fastq_it0(line_iterator0),
                fastq_it1(line_iterator1);

  size_t seqnum = 0, total_length[2] = {0};
  auto it0 = fastq_it0.begin();
  auto it1 = fastq_it1.begin();
  while (it0 != fastq_it0.end() && it1 != fastq_it1.end())
  {
    const std::string_view &sequence0 = (*it0).sequence_get();
    const std::string_view &sequence1 = (*it1).sequence_get();
    if (fasta_output)
    {
      const std::string_view &header0 = (*it0).header_get();
      const std::string_view &header1 = (*it1).header_get();
      std::cout << ">" << header0.substr(1) << std::endl;
      std::cout << sequence0 << std::endl;
      std::cout << ">" << header1.substr(1) << std::endl;
      std::cout << sequence1 << std::endl;
    }
    total_length[0] += sequence0.size();
    total_length[1] += sequence1.size();
    seqnum++;
    ++it0;
    ++it1;
  }
  if (statistics)
  {
    std::cout << "# number of sequences\t" << seqnum << std::endl;
    std::cout << "# total length0\t" << total_length[0] << std::endl;
    std::cout << "# total length1\t" << total_length[1] << std::endl;
  }
}

static void parallel_char_distribution(
      const char *file_contents,
      const std::vector<std::pair<size_t,size_t>> &intervals)
{
  size_t *count_entries = static_cast<size_t *>(calloc(intervals.size(),
                                                       sizeof *count_entries)),
         *dist = static_cast<size_t *>(calloc(4 * intervals.size(),
                                              sizeof *dist));
  std::vector<std::thread> threads{};
  size_t thd_num = 0;
  for (auto &&itv : intervals)
  {
    threads.push_back(std::thread([file_contents,&itv,count_entries,dist,
                                   thd_num]
    {
      GttlLineIterator<0> line_iterator(file_contents + std::get<0>(itv),
                                        std::get<1>(itv));
      GttlFastQIterator<GttlLineIterator<0>> fastq_it(line_iterator);
      size_t local_count_entries = 0;
      size_t *local_dist = dist + 4 * thd_num;
      for (auto &&fastq_entry : fastq_it)
      {
        local_count_entries++;
        const std::string_view &sequence = fastq_entry.sequence_get();
        for (auto &&cc : sequence)
        {
          local_dist[(static_cast<uint8_t>(cc) >> 1) & uint8_t(3)]++;
        }
      }
      count_entries[thd_num] = local_count_entries;
    }));
    thd_num++;
  }
  for (auto &th : threads)
  {
    th.join();
  }
  size_t total_count_entries = 0;
  for (size_t thd_num = 0; thd_num < intervals.size(); thd_num++)
  {
    total_count_entries += count_entries[thd_num];
  }
  std::cout << "# total_count_entries\t" << total_count_entries << std::endl;
  free(count_entries);
  for (int char_idx = 0; char_idx < 4; char_idx++)
  {
    size_t total_cc_count = 0;
    for (size_t thd_num = 0; thd_num < intervals.size(); thd_num++)
    {
      total_cc_count += dist[4 * thd_num + char_idx];
    }
    std::cout << "# char\t" << char_idx << "\t" << total_cc_count << std::endl;
  }
  free(dist);
}

int main(int argc,char *argv[])
{
  SeqReaderOptions options{2,true};

  try
  {
    options.parse(argc, argv);
  }
  catch (std::invalid_argument &e) /* check_err.py */
  {
    std::cerr << argv[0] << ": " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  if (options.help_option_is_set())
  {
    return EXIT_SUCCESS;
  }
  try
  {
    const bool statistics = options.statistics_option_is_set();
    const bool echo = options.echo_option_is_set();
    const bool fasta_output = options.fasta_output_option_is_set();
    const size_t split_size = options.split_size_get();
    const std::vector<std::string> &inputfiles = options.inputfiles_get();
    if (inputfiles.size() == 1)
    {
      if (options.num_threads_get() > 0)
      {
        FastQParts fastq_parts(options.num_threads_get(),inputfiles[0]);
        fastq_parts.show();
        parallel_char_distribution(fastq_parts.file_contents,
                                   fastq_parts.intervals);
      } else
      {
        if (split_size > 0)
        {
          fastq_split_writer(split_size,inputfiles[0]);
        } else
        {
          if (options.mapped_option_is_set())
          {
            process_single_file_mapped(statistics,echo,fasta_output,
                                       inputfiles[0]);
          } else
          {
            process_single_file_streamed(statistics,echo,fasta_output,
                                         inputfiles[0]);
          }
        }
      }
    } else
    {
      process_paired_files(statistics,fasta_output,inputfiles[0],inputfiles[1]);
    }
  }
  catch (std::string &msg)
  {
    std::cerr << argv[0] << ": " << msg << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
