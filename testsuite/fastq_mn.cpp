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
#include <iomanip>
#include <cassert>
#include <algorithm>
#include <thread>
#include <limits>
#include <map>
#include "utilities/mathsupport.hpp"
#include "utilities/str_format.hpp"
#include "utilities/basename.hpp"
#include "utilities/gttl_mmap.hpp"
#include "utilities/wyhash.hpp"
#include "utilities/mathsupport.hpp"
#include "utilities/file_size.hpp"
#include "threading/threadsafe_queue.hpp"
#ifdef WITH_XXHASH
#define XXH_INLINE_ALL
#include "utilities/xxhash.hpp"
#endif
#include "sequences/gttl_fastq_generator.hpp"
#include "sequences/split.hpp"
#include "sequences/dna_seq_encoder.hpp"
#include "sequences/dna_seq_decoder.hpp"
#include "sequences/format_sequence.hpp"
#include "seq_reader_options.hpp"

static void fastq_split_writer(size_t split_size,
                               const std::string &inputfilename)
{
  constexpr const int buf_size = 1 << 14;
  GttlFastQGenerator<buf_size> fastq_it(inputfilename.c_str());
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
      const std::string_view &sequence = (*it)->sequence_get();
      const std::string_view &header = (*it)->header_get();
      const std::string_view &quality = (*it)->quality_get();
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
                               size_t fasta_output,
                               hash_mode_type hash_mode,
                               const std::string &inputfilename,
                               FastQIterator &fastq_it)
{
  size_t seqnum = 0,
         total_length = 0,
         min_length = std::numeric_limits<size_t>::max(),
         max_length = 0;;
  uint64_t hash_value_sum = 0;
  std::map<size_t,size_t> length_dist_map{};
  for (auto &&fastq_entry : fastq_it)
  {
    const std::string_view &sequence = fastq_entry->sequence_get();
    if (statistics)
    {
      min_length = std::min(min_length,sequence.size());
      max_length = std::max(max_length,sequence.size());
      length_dist_map[sequence.size()]++;
    }
    if (hash_mode == hash_mode_wy)
    {
      const uint64_t hash_value = wyhash(sequence.data(), sequence.size(),0);
      if (echo)
      {
        std::cout << hash_value << "\t" << sequence << std::endl;
      } else
      {
        hash_value_sum += hash_value;
      }
    } else
    {
      if (hash_mode == hash_mode_xx)
      {
#ifdef WITH_XXHASH
        const uint64_t hash_value = XXH64(sequence.data(),sequence.size(),0);
#else
        const uint64_t hash_value = 0;
#endif
        if (echo)
        {
          std::cout << hash_value << "\t" << sequence << std::endl;
        } else
        {
          hash_value_sum += hash_value;
        }
      } else
      {
        if (echo)
        {
          const std::string_view &header = fastq_entry->header_get();
          const std::string_view &quality = fastq_entry->quality_get();
          std::cout << header << std::endl
                    << sequence << std::endl
                    << "+" << std::endl
                    << quality << std::endl;
        } else
        {
          if (fasta_output > 0)
          {
            const std::string_view &header = fastq_entry->header_get();
            std::cout << ">" << header.substr(1) << std::endl;
            gttl_format_sequence(sequence,fasta_output);
          }
        }
      }
    }
    total_length += sequence.size();
    seqnum++;
  }
  if (statistics)
  {
    std::cout << "# number of sequences\t" << seqnum << std::endl;
    std::cout << "# total length\t" << total_length << std::endl;
    std::cout << "# mean sequence length\t" << total_length/seqnum << std::endl;
    std::cout << "# minimum sequence length\t" << min_length << std::endl;
    std::cout << "# maximum sequence length\t" << max_length << std::endl;
    std::cout << "# length, count, count ratio" << std::endl;
    for (auto const& [len_as_key, count_as_value] : length_dist_map)
    {
      std::cout << "# " << len_as_key << "\t" << count_as_value << "\t"
                << std::fixed << std::setprecision(4)
                << (100.0 * static_cast<double>(count_as_value)/seqnum)
                << std::endl;
    }
    std::cout << "# original size of file " << inputfilename << " (MB)\t"
              << static_cast<size_t>(mega_bytes(gttl_file_size(inputfilename)))
              << std::endl;
    std::cout << "# expected size of sequences in RAM (MB)\t"
              << static_cast<size_t>(mega_bytes(total_length / 4))
              << std::endl;
    const size_t uint64_t_units = ((max_length * 2) + 63)/64;
    std::cout << "# expected size of DNAEncoding in RAM (MB)\t"
              << static_cast<size_t>(mega_bytes(seqnum * 8 * uint64_t_units))
              << std::endl;
  }
  if (hash_mode != hash_mode_none)
  {
    std::cout << inputfilename << "\t" << hash_value_sum << std::endl;
  }
}

static void process_single_file_streamed(bool statistics,
                                         bool echo,
                                         size_t fasta_output,
                                         hash_mode_type hash_mode,
                                         const std::string &inputfilename)
{
  constexpr const int buf_size = 1 << 14;
  using FastQIterator = GttlFastQGenerator<buf_size>;
  FastQIterator fastq_it(inputfilename.c_str());
  process_fastq_iter<FastQIterator>(statistics,echo,fasta_output,hash_mode,
                                    inputfilename,fastq_it);
}

static void process_single_file_mapped(bool statistics,
                                       bool echo,
                                       size_t fasta_output,
                                       hash_mode_type hash_mode,
                                       const std::string &inputfilename)
{
  constexpr const int buf_size = 1 << 14;
  Gttlmmap<char> mapped_file(inputfilename.c_str());
  using FastQIterator = GttlFastQGenerator<buf_size>;
  FastQIterator fastq_it(mapped_file.ptr(), mapped_file.size());
  process_fastq_iter<FastQIterator>(statistics,echo,fasta_output,hash_mode,
                                    inputfilename,fastq_it);
}

static void process_paired_files(bool statistics,
                                 size_t fasta_output,
                                 const std::string &filename0,
                                 const std::string &filename1)
{
  constexpr const int buf_size = 1 << 14;
  using FastQIterator = GttlFastQGenerator<buf_size>;
  FastQIterator fastq_it0(filename0.c_str()),
                fastq_it1(filename1.c_str());

  size_t seqnum = 0, total_length[2] = {0};
  auto it0 = fastq_it0.begin();
  auto it1 = fastq_it1.begin();
  while (it0 != fastq_it0.end() && it1 != fastq_it1.end())
  {
    const std::string_view &sequence0 = (*it0)->sequence_get();
    const std::string_view &sequence1 = (*it1)->sequence_get();
    if (fasta_output > 0)
    {
      const std::string_view &header0 = (*it0)->header_get();
      const std::string_view &header1 = (*it1)->header_get();
      std::cout << ">" << header0.substr(1) << std::endl;
      gttl_format_sequence(sequence0,fasta_output);
      std::cout << ">" << header1.substr(1) << std::endl;
      gttl_format_sequence(sequence1,fasta_output);
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

static void char_distribution_seq(const std::string &inputfilename)
{
  GttlFpType in_fp = gttl_fp_type_open(inputfilename.c_str(), "rb");
  if (in_fp == nullptr)
  {
    throw std::string(": cannot open file");
    /* check_err.py checked */
  }
  constexpr const int buf_size = 1 << 14;
  GttlFastQGenerator<buf_size> fastq_it(in_fp);
  size_t dist[4] = {0};
  size_t count_entries = 0;
  for (auto fastq_entry : fastq_it)
  {
    count_entries++;
    const std::string_view &sequence = fastq_entry->sequence_get();
    for (auto &&cc : sequence)
    {
      dist[(static_cast<uint8_t>(cc) >> 1) & uint8_t(3)]++;
    }
  }
  gttl_fp_type_close(in_fp);
  std::cout << "# total_count_entries\t" << count_entries << std::endl;
  for (int char_idx = 0; char_idx < 4; char_idx++)
  {
    std::cout << "# char\t" << char_idx << "\t" << dist[char_idx] << std::endl;
  }
}

static void char_distribution_thd_gz(size_t num_threads,
                                     const std::string &inputfilename)
{
  assert(num_threads >= 2);
  GttlFpType in_fp = gttl_fp_type_open(inputfilename.c_str(), "rb");
  if (in_fp == nullptr)
  {
    throw std::string(": cannot open file");
    /* check_err.py checked */
  }
  constexpr const int buf_size = 1 << 12;
  using BufferedFastQIter = GttlFastQGenerator<buf_size>;
  BufferedFastQIter fastq_it(in_fp);
  size_t *dist = static_cast<size_t *>(calloc(4 * (num_threads-1),
                                              sizeof *dist));
  std::vector<std::thread> threads{};
  ThreadsafeQueue<std::string> sequence_queue;
  threads.push_back(std::thread([&sequence_queue, &fastq_it]
  {
    size_t count_entries = 0;
    for (auto &&fastq_entry : fastq_it)
    {
      count_entries++;
      const std::string_view &seq_view = fastq_entry->sequence_get();
      sequence_queue.enqueue(std::string(seq_view.begin(),seq_view.end()));
    }
    std::cout << "# total_count_entries\t" << count_entries << std::endl;
  }));
  for (size_t thd_num = 1; thd_num < num_threads; thd_num++)
  {
    threads.push_back(std::thread([&sequence_queue, &dist, thd_num]
    {
      size_t *local_dist = dist + 4 * (thd_num-1);
      using namespace std::chrono_literals;
      std::this_thread::sleep_for(400ms);
      while (true)
      {
        std::optional<std::string> opt_sequence = sequence_queue.dequeue();
        if (not opt_sequence)
        {
          break;
        }
        auto sequence = opt_sequence.value();
        for (auto cc : sequence)
        {
          local_dist[(static_cast<uint8_t>(cc) >> 1) & uint8_t(3)]++;
        }
      }
    }));
  }
  for (auto &th : threads)
  {
    th.join();
  }
  gttl_fp_type_close(in_fp);
  size_t flushed_sequences = 0;
  while (true)
  {
    std::optional<std::string> opt_sequence = sequence_queue.dequeue();
    if (not opt_sequence)
    {
      break;
    }
    auto sequence = opt_sequence.value();
    for (auto cc : sequence)
    {
      dist[(static_cast<uint8_t>(cc) >> 1) & uint8_t(3)]++;
    }
    flushed_sequences++;
  }
  std::cout << "# flushed_sequences\t" << flushed_sequences << std::endl;
  for (int char_idx = 0; char_idx < 4; char_idx++)
  {
    size_t total_cc_count = 0;
    for (size_t thd_num = 1; thd_num < num_threads; thd_num++)
    {
      total_cc_count += dist[4 * (thd_num-1) + char_idx];
    }
    std::cout << "# char\t" << char_idx << "\t" << total_cc_count << std::endl;
  }
  free(dist);
}

static void char_distribution_thd(const SequencesSplit &sequences_split)
{
  size_t *count_entries = static_cast<size_t *>(calloc(sequences_split.size(),
                                                       sizeof *count_entries)),
         *dist = static_cast<size_t *>(calloc(4 * sequences_split.size(),
                                              sizeof *dist));
  std::vector<std::thread> threads{};
  for (size_t thd_num = 0; thd_num < sequences_split.size(); thd_num++)
  {
    threads.push_back(std::thread([&sequences_split, count_entries,dist,thd_num]
    {
      const std::string_view &this_view = sequences_split[thd_num];
      GttlFastQGenerator<16384> fastq_it(this_view.data(), this_view.size());
      size_t local_count_entries = 0;
      size_t *local_dist = dist + 4 * thd_num;
      for (auto &&fastq_entry : fastq_it)
      {
        local_count_entries++;
        const std::string_view &sequence = fastq_entry->sequence_get();
        for (auto &&cc : sequence)
        {
          local_dist[(static_cast<uint8_t>(cc) >> 1) & uint8_t(3)]++;
        }
      }
      count_entries[thd_num] = local_count_entries;
    }));
  }
  for (auto &th : threads)
  {
    th.join();
  }
  size_t total_count_entries = 0;
  for (size_t thd_num = 0; thd_num < sequences_split.size(); thd_num++)
  {
    total_count_entries += count_entries[thd_num];
  }
  std::cout << "# total_count_entries\t" << total_count_entries << std::endl;
  free(count_entries);
  for (int char_idx = 0; char_idx < 4; char_idx++)
  {
    size_t total_cc_count = 0;
    for (size_t thd_num = 0; thd_num < sequences_split.size(); thd_num++)
    {
      total_cc_count += dist[4 * thd_num + char_idx];
    }
    std::cout << "# char\t" << char_idx << "\t" << total_cc_count << std::endl;
  }
  free(dist);
}

static std::string qgram_decode(uint64_t code,size_t qgram_length)
{
  static const std::array<char,4> dna_letters{'A','C','G','T'};
  std::string s(qgram_length,' ');
  for (int idx = qgram_length - 1; idx >= 0; idx--)
  {
    s[idx] = dna_letters[code & uint64_t(3)];
    code >>= 2;
  }
  return s;
}

static void verify_consecutive_qgrams(const uint64_t *sub_unit_ptr,
                                      size_t qgram_length,
                                      size_t sequence_length)
{
  assert(qgram_length <= sequence_length);
  DNAQgramDecoder
    dna_qgram_decoder(qgram_length,
                      sequence_length + 1 - qgram_length,
                      sub_unit_ptr);
  std::string previous_qgram;
  for (auto const qgram_code : dna_qgram_decoder)
  {
    const std::string qgram(qgram_decode(qgram_code,qgram_length));
    //std::cout << qgram_code << "\t"
              //<< qgram << std::endl;
    if (previous_qgram.size() > 0)
    {
      for (size_t idx = 0; idx < qgram_length-1; idx++)
      {
        const char p_cc = previous_qgram[idx+1],
                     cc = qgram[idx];
        if (p_cc != cc)
        {
          std::cerr << "p_cc = " << p_cc << " != " << cc
                    << std::endl;
          exit(EXIT_FAILURE);
        }
      }
    }
    previous_qgram = qgram;
  }
}

template<bool split_at_wildcard>
static void verify_decoding_multilength(bool statistics,
                                        const std::string &inputfilename,
                                        size_t qgram_length)
{
  DNAEncodingMultiLength<uint64_t,split_at_wildcard,false>
    dna_encoding_multi_length(inputfilename);
  if (statistics)
  {
    dna_encoding_multi_length.statistics();
  }
  dna_encoding_multi_length.prepare_view(1);
  size_t seqcount = 0;
  for (auto const &[sub_unit_ptr, sequence_length] : dna_encoding_multi_length)
  {
    assert(sub_unit_ptr != nullptr);
    verify_consecutive_qgrams(sub_unit_ptr,qgram_length,sequence_length);
    seqcount++;
  }
  std::cout << "# verified " << qgram_length << "-mers in " << seqcount
            << " sequences" << std::endl;
}

static void verify_decoding_parts_view(
  const DNAEncodingMultiLength<uint64_t,false,false> &dna_encoding_multi_length)
{
  std::map<size_t,size_t> length_dist_map;
  for (size_t part_idx = 0;
       part_idx < dna_encoding_multi_length.num_parts_get();
       part_idx++)
  {
    for(auto it = dna_encoding_multi_length.begin(part_idx);
        it != dna_encoding_multi_length.end(); ++it)
    {
      length_dist_map[std::get<1>(*it)]++;
    }
  }
  dna_encoding_multi_length.verify_length_dist(length_dist_map);
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
    const size_t fasta_output = options.fasta_output_option_is_set()
                                  ? options.line_width_get() : 0;
    const hash_mode_type hash_mode = options.hash_mode_get();
    const size_t split_size = options.split_size_get();
    const std::vector<std::string> &inputfiles = options.inputfiles_get();
    if (inputfiles.size() == 1)
    {
      if (options.num_threads_get() > 0)
      {
        if (options.num_threads_get() == 1)
        {
          char_distribution_seq(inputfiles[0]);
        } else
        {
          if (gttl_has_suffix(inputfiles[0],std::string(".gz")))
          {
            try
            {
              char_distribution_thd_gz(options.num_threads_get(),
                                       inputfiles[0]);
            }
            catch (std::string &msg)
            {
              std::cerr << argv[0] << ": " << msg << std::endl;
              return EXIT_FAILURE;
            }
          } else
          {
            const bool fasta_format = false;
            SequencesSplit sequences_split(options.num_threads_get(),
                                           inputfiles[0],
                                           fasta_format);
            sequences_split.show();
            char_distribution_thd(sequences_split);
          }
        }
      } else
      {
        if (split_size > 0)
        {
          fastq_split_writer(split_size,inputfiles[0]);
        } else
        {
          static constexpr bool split_at_wildcard = false;
          if (options.encoding_type_get() != std::string(""))
          {
            if (options.encoding_type_get() == std::string("8"))
            {
              DNAEncodingMultiLength<uint8_t,split_at_wildcard,false>
                dna_encoding_multi_length(inputfiles[0]);
              if (statistics)
              {
                dna_encoding_multi_length.statistics();
              }
            } else
            {
              if (options.encoding_type_get() == std::string("16"))
              {
                DNAEncodingMultiLength<uint16_t,split_at_wildcard,false>
                  dna_encoding_multi_length(inputfiles[0]);
                if (statistics)
                {
                  dna_encoding_multi_length.statistics();
                }
              } else
              {
                if (options.encoding_type_get() == std::string("32"))
                {
                  DNAEncodingMultiLength<uint32_t,split_at_wildcard,false>
                    dna_encoding_multi_length(inputfiles[0]);
                  if (statistics)
                  {
                    dna_encoding_multi_length.statistics();
                  }
                } else
                {
                  if (options.encoding_type_get() == std::string("64"))
                  {
                    DNAEncodingMultiLength<uint64_t,split_at_wildcard,false>
                      dna_encoding_multi_length(inputfiles[0]);
                    if (statistics)
                    {
                      dna_encoding_multi_length.statistics();
                    } else
                    {
                      for (size_t num_parts = 1;
                           num_parts <
                              std::min(dna_encoding_multi_length.
                                       total_number_of_nucleotide_ranges_get(),
                                       size_t(10));
                           num_parts++)
                      {
                        dna_encoding_multi_length.prepare_view(num_parts);
                        verify_decoding_parts_view
                           (dna_encoding_multi_length);
                      }
                    }
                  } else
                  {
                    int bits, r_qgram_length;
                    if (std::sscanf(options.encoding_type_get().c_str(),"%d,%d",
                                    &bits,&r_qgram_length) != 2
                        or bits != 64
                        or r_qgram_length < 2
                        or r_qgram_length > 32)
                    {
                      std::cerr << argv[0]
                                << ": argument of --encoding must be one of"
                                << " 8, 16, 32, 64,<qgram_length>"
                                << std::endl;
                      return EXIT_FAILURE;
                    }
                    verify_decoding_multilength<split_at_wildcard>
                                               (statistics,
                                                inputfiles[0],
                                                static_cast<size_t>
                                                           (r_qgram_length));
                  }
                }
              }
            }
          } else
          {
            if (options.mapped_option_is_set())
            {
              process_single_file_mapped(statistics,echo,fasta_output,
                                         hash_mode,inputfiles[0]);
            } else
            {
              process_single_file_streamed(statistics,echo,fasta_output,
                                           hash_mode,inputfiles[0]);
            }
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
    const std::vector<std::string> &inputfiles = options.inputfiles_get();
    for (auto &&inputfile : inputfiles)
    {
      std::cerr << argv[0] << ": file \"" << inputfile << "\""
                << msg << std::endl;
    }
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
