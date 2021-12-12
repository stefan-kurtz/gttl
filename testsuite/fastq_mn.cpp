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
#include <cstdbool>
#include <cassert>
#include <algorithm>
#include "sequences/gttl_fastq_iterator.hpp"

typedef enum
{
  EchoComplete,
  EchoRead,
  Statistics
} OutputMode;

struct SeqPairCount
{
  size_t num_sequences = 0,
         total_length[2] = {0};
};

void call_back_seq_pair_count(SeqPairCount *spc,
                              const std::string *seq1,
                              const std::string *seq2)
{
  spc->num_sequences++;
  spc->total_length[0] += seq1->size();
  spc->total_length[1] += seq2->size();
}

static std::pair<size_t,size_t> process_single_file(OutputMode output_mode,
                                                    const char *filename)
{
  constexpr const int buf_size = 1 << 14;
  GttlFastQIterator<buf_size> fastq_it(filename);
  size_t seqnum = 0, total_length = 0;

  for (auto &&fastq_entry : fastq_it)
  {
    const std::string &sequence = std::get<1>(fastq_entry);
    total_length += sequence.size();
    if (output_mode == EchoComplete)
    {
      const std::string &quality = std::get<3>(fastq_entry);
      std::cout << fastq_header(fastq_entry) << std::endl;
      std::cout << sequence << std::endl;
      std::cout << "+" << std::endl;
      std::cout << quality << std::endl;
    } else
    {
      if (output_mode == EchoRead)
      {
        std::cout << ">" << seqnum << std::endl;
        std::cout << sequence << std::endl;
      }
    }
    seqnum++;
  }
  return {seqnum, total_length};
}

int main(int argc,char *argv[])
{
  if ((argc != 3 && argc != 4) ||
      (strcmp(argv[1],"E") != 0 &&
       strcmp(argv[1],"S") != 0 &&
       strcmp(argv[1],"R") != 0) ||
      (argc == 4 && (strcmp(argv[1],"E") == 0 || strcmp(argv[1],"R") == 0)))
  {
    std::cerr << "Usage: " << argv[0] << " E|S|R <inputfile> [inputfile]"
              << std::endl;
    std::cerr << "  S: output statistics"
              << std::endl;
    std::cerr << "  E: echo input (only available when processing single file)"
              << std::endl;
    std::cerr << "  R: echo only read (only available when processing "
                 "single file)"
              << std::endl;
    return EXIT_FAILURE;
  }
  OutputMode output_mode = EchoComplete;
  if (strcmp(argv[1],"S") == 0)
  {
    output_mode = Statistics;
  } else
  {
    if (strcmp(argv[1],"E") == 0)
    {
      output_mode = EchoComplete;
    } else
    {
      output_mode = EchoRead;
    }
  }
  const char *progname = argv[0];
  const char *filename1 = argv[2];
  std::pair<size_t,size_t> result;
  SeqPairCount seq_pair_count{};
  try
  {
    if (argc == 3)
    {
      result = process_single_file(output_mode,filename1);
    } else
    {
      assert(argc == 4);
#ifdef READPAIRS
      const char *filename2 = argv[3];
      quick_seq_pair_reader<SeqPairCount,call_back_seq_pair_count>
                           (&seq_pair_count,filename1,filename2);
#endif
    }
  }
  catch (std::string &msg)
  {
    std::cerr << progname << ": " << msg << std::endl;
    return EXIT_FAILURE;
  }
  if (output_mode == Statistics)
  {
    if (argc == 3)
    {
      size_t num_sequences = std::get<0>(result);
      size_t total_length = std::get<1>(result);
      std::cout << "# number of sequences\t" << num_sequences << std::endl;
      std::cout << "# total length\t" << total_length << std::endl;
      std::cout << "# mean length\t" << total_length/num_sequences << std::endl;
    } else
    {
      assert(argc == 4);
      std::cout << "# number of sequences\t" << seq_pair_count.num_sequences
                << std::endl;
      std::cout << "# total length1\t" << seq_pair_count.total_length[0]
                << std::endl;
      std::cout << "# mean length\t" << seq_pair_count.total_length[0]/
                                        seq_pair_count.num_sequences
                << std::endl;
    }
  }
  return EXIT_SUCCESS;
}
