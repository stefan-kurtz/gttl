#ifndef SPLIT_FILES_HPP
#define SPLIT_FILES_HPP

#include <cstddef>
#include <cmath>
#include <functional>
#include <sstream>
#include <string>
#include <string_view>
#include "threading/thread_pool_unknown_tasks.hpp"
#include "utilities/write_output_file.hpp"

/*
** Split a FastQGenerator or FastAGenerator into fragments of a given length (of
** symbols).
** base_name describes the base of the output filename. File extensions and an
** increasing
** number will be added onto this.
** Calling this function with `part_length=0` will produce single-sequence
** files.
** This is because any non-empty sequence will always be longer than 0
** characters.
*/
template <class SequenceGeneratorClass>
void split_into_parts_length(SequenceGeneratorClass &seq_gen,
                             const std::string &base_name,
                             size_t part_length,
                             size_t compression_level,
                             size_t n_threads,
                             size_t padding_length = 2)
{
  size_t part_number = 1;
  size_t length_iterated = 0;
  std::ostringstream s_out;
  ThreadPoolUnknownTasks<void, void> tp(n_threads);
  const std::string output_file_suffix{SequenceGeneratorClass::
                                         is_fastq_generator ? ".fastq"
                                                            : ".fasta"};

  for (const auto *si : seq_gen)
  {
    const std::string_view &sequence = si->sequence_get();
    const std::string_view &header = si->header_get();

    if constexpr (SequenceGeneratorClass::is_fastq_generator)
    {
      s_out << "@";
    } else
    {
      s_out << ">";
    }
    s_out << header << '\n' << sequence << '\n';

    if constexpr (SequenceGeneratorClass::is_fastq_generator)
    {
      const std::string_view &quality = si->quality_get();
      s_out << "+\n" << quality << '\n';
    }

    length_iterated += sequence.size();

    if (length_iterated >= part_length)
    {
      std::string fname_out = base_name;
      for(size_t i = 1; i <= padding_length; i++)
      {
        fname_out += (part_number < static_cast<size_t>(std::pow(10, i))
                        ? "0"
                        : "");
      }
      fname_out += std::to_string(part_number) + output_file_suffix;
      if (n_threads == 1)
      {
        write_to_output_file(fname_out, s_out.str(), compression_level);
      } else
      {
        tp.enqueue([fname_out, capture0 = s_out.str(), compression_level] {
          write_to_output_file(fname_out, capture0, compression_level);
        });
      }

      s_out.str("");
      length_iterated = 0;
      part_number++;
    }
  }
  if (not s_out.str().empty())
  {
    std::string fname_out = base_name;
    for(size_t i = 1; i <= padding_length; i++)
    {
      fname_out += (part_number < static_cast<size_t>(std::pow(10, i))
                    ? "0"
                    : "");
    }
    fname_out += std::to_string(part_number) + output_file_suffix;
    if (n_threads == 1)
    {
      write_to_output_file(fname_out, s_out.str(), compression_level);
    } else
    {
      tp.enqueue([fname_out, capture0 = s_out.str(), compression_level] {
        write_to_output_file(fname_out, capture0, compression_level);
      });
    }
  }
}

/*
** Split a FastQGenerator or FastAGenerator into fragments of a given number of
** sequences each.
*/
template <class SequenceGeneratorClass>
void split_into_num_sequences(SequenceGeneratorClass &seq_gen,
                              const std::string &base_name,
                              size_t seqs_per_file,
                              size_t compression_level,
                              size_t n_threads)
{
  ThreadPoolUnknownTasks<void,void> tp(n_threads);
  size_t part_number = 1;
  size_t seqs_iterated = 0;
  const std::string output_file_suffix{SequenceGeneratorClass::
                                         is_fastq_generator ? ".fastq"
                                                            : ".fasta"};
  std::ostringstream s_out;
  for (const auto *si : seq_gen)
  {
    const std::string_view &sequence = si->sequence_get();
    const std::string_view &header = si->header_get();

    if constexpr (SequenceGeneratorClass::is_fastq_generator)
    {
      s_out << "@";
    } else
    {
      s_out << ">";
    }

    s_out << header << '\n' << sequence << '\n';

    if constexpr (SequenceGeneratorClass::is_fastq_generator)
    {
      const std::string_view &quality = si->quality_get();
      s_out << "+\n" << quality << '\n';
    }

    seqs_iterated++;

    if (seqs_iterated >= seqs_per_file)
    {
      const std::string fname_out = base_name + (part_number <= 9 ? "0" : "")
                                   + std::to_string(part_number)
                                   + output_file_suffix;
      if (n_threads == 1)
      {
        write_to_output_file(fname_out, s_out.str(), compression_level);
      } else
      {
        tp.enqueue([fname_out, capture0 = s_out.str(), compression_level] {
          write_to_output_file(fname_out, capture0, compression_level);
        });
      }

      s_out.str("");
      seqs_iterated = 0;
      part_number++;
    }
  }
  if (not s_out.str().empty())
  {
    const std::string fname_out = base_name + (part_number <= 9 ? "0" : "")
                                 + std::to_string(part_number)
                                 + output_file_suffix;
    if (n_threads == 1)
    {
      write_to_output_file(fname_out, s_out.str(), compression_level);
    } else
    {
      tp.enqueue([fname_out, capture0 = s_out.str(), compression_level] {
        write_to_output_file(fname_out, capture0, compression_level);
      });
    }
  }
}

/*
** Split a FastQGenerator or FastAGenerator into a given number of fragments.
** This is a wrapper around split_into_parts_length which will iterate
** over the entire input file once and determine total sequence length and thus
** the necessary length of each sequence part.
*/
template <class SequenceGeneratorClass>
void split_into_num_files(SequenceGeneratorClass &seq_gen,
                          const std::string &base_name,
                          size_t part_num,
                          size_t compression_level,
                          size_t n_threads)
{
  size_t total_length = 0;
  for (const auto *si : seq_gen)
  {
    total_length += si->sequence_get().size();
  }
  assert(part_num > 0);
  const size_t part_len = (total_length + part_num - 1)/ part_num;
  seq_gen.reset();
  split_into_parts_length(seq_gen, base_name, part_len, compression_level,
                          n_threads, static_cast<size_t>(std::log10(part_num)));
}

#endif // SPLIT_FILES_HPP
