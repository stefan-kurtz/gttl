#ifndef SPLIT_FILES_HPP
#define SPLIT_FILES_HPP

#include <cstddef>
#include <sstream>
#include <type_traits>
#include "sequences/gttl_fastq_iterator.hpp"
#include "utilities/write_output_file.hpp"

template <class T>
struct is_fastq_iterator : std::false_type {};

template <class T>
struct is_fastq_iterator<GttlFastQIterator<T>> : std::true_type {};


/*
** Split a FastQIterator or SeqIterator into fragments of a given length (of
** symbols).
** base_name describes the base of the output filename. File extensions and an
** increasing
** number will be added onto this.
** Calling this function with `part_length=0` will produce single-sequence
** files.
** This is because any non-empty sequence will always be longer than 0
** characters.
*/
template <class SequenceIterator>
void split_into_parts_length(SequenceIterator &seq_it,
                             const std::string &base_name,
                             const size_t &part_length,
                             const size_t &compression_level)
{
  size_t part_number = 1, length_iterated = 0;
  std::ostringstream s_out;
  for (auto &&si : seq_it)
  {
    const std::string_view &sequence = si.sequence_get();
    const std::string_view &header = si.header_get();

    s_out << header << "\n" << sequence << "\n";

    if constexpr (is_fastq_iterator<SequenceIterator>::value)
    {
      const std::string_view &quality = si.quality_get();
      s_out << "+\n" << quality << "\n";
    }

    length_iterated += sequence.size();

    if (length_iterated > part_length)
    {
      std::string fname_out = base_name + std::to_string(part_number) +
                              (seq_it.is_fastq_iterator ? ".fastq" : ".fasta");
      write_to_output_file(fname_out, s_out.str(), compression_level);

      s_out.str("");
      length_iterated = 0;
      part_number++;
    }
  }
  if (!s_out.str().empty())
  {
    std::string fname_out = base_name + std::to_string(part_number) +
                            (seq_it.is_fastq_iterator ? ".fastq" : ".fasta");
    write_to_output_file(fname_out, s_out.str(), compression_level);
  }
}


/*
** Split a FastQIterator or SeqIterator into fragments of a given number of
** sequences each.
*/
template <class SequenceIterator>
void split_into_num_sequences(SequenceIterator &seq_it,
                              const std::string &base_name,
                              size_t seqs_per_file, size_t compression_level)
{
  size_t part_number = 1, seqs_iterated = 0;
  std::ostringstream s_out;
  for (auto &&si : seq_it)
  {
    const std::string_view &sequence = si.sequence_get();
    const std::string_view &header = si.header_get();

    s_out << header << "\n" << sequence << "\n";

    if constexpr (seq_it.is_fastq_iterator)
    {
      const std::string_view &quality = si.quality_get();
      s_out << "+\n" << quality << "\n";
    }

    seqs_iterated++;

    if (seqs_iterated >= seqs_per_file)
    {
      std::string fname_out = base_name + std::to_string(part_number) +
                              (seq_it.is_fastq_iterator ? ".fastq" : ".fasta");
      write_to_output_file(fname_out, s_out.str(), compression_level);

      s_out.str("");
      seqs_iterated = 0;
      part_number++;
    }
  }
  if (!s_out.str().empty())
  {
    std::string fname_out = base_name + std::to_string(part_number) +
                            (seq_it.is_fastq_iterator ? ".fastq" : ".fasta");
    write_to_output_file(fname_out, s_out.str(), compression_level);
  }
}


/*
** Split a FastQIterator or SeqIterator into a given number of fragments.
** This is a wrapper around split_into_parts_length which will iterate
** over the entire input file once and determine total sequence length and thus
** the necessary length of each sequence part.
*/
template <class SequenceIterator>
void split_into_num_files(SequenceIterator &seq_it,
                          const std::string &base_name, size_t part_num,
                          size_t compression_level)
{
  size_t total_length = 0;
  for (auto &&si : seq_it)
  {
    total_length += si.sequence_get().size();
  }
  const size_t part_len =
      total_length / part_num + (size_t) (total_length % part_num != 0);
  seq_it.reset();
  split_into_parts_length(seq_it, base_name, part_len, compression_level);
}

#endif // SPLIT_FILES_HPP
