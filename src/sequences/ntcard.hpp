#ifndef NTCARD_HPP
#define NTCARD_HPP

#include <thread>
#include <vector>
#include "utilities/gttl_file_open.hpp"
#include "utilities/has_suffix_or_prefix.hpp"
#include "sequences/char_range.hpp"
#include "sequences/char_finder.hpp"
#include "sequences/qgrams_hash_nthash.hpp"
#include "sequences/gttl_seq_iterator.hpp"
#include "sequences/gttl_fastq_iterator.hpp"
#include "sequences/split.hpp"
#include "sequences/dna_seq_encoder.hpp"
#include "sequences/dna_seq_decoder.hpp"

template <bool split_at_wildcard,
          class SeqIterator,
          class HashValueIterator,
          class TableClass>
static void ntcard_enumerate_inner(SeqIterator* gttl_si,
                                   TableClass* table,
                                   size_t qgram_length)
{
  using NucleotideRanger = GttlCharRange<char_finder::NucleotideFinder,
                                         ntcard_nucleotide_finder, true, false>;

  for (auto &&si : *gttl_si)
  {
    auto sequence = si.sequence_get();
    if constexpr (split_at_wildcard)
    {
      NucleotideRanger nuc_ranger(sequence.data(), sequence.size());
      for (auto const &&range : nuc_ranger)
      {
        const size_t this_length = std::get<1>(range);
        if (this_length >= qgram_length)
        {
          const char *substring = sequence.data() + std::get<0>(range);
          HashValueIterator qgiter(qgram_length, substring, this_length);
          for (auto const &&code_pair : qgiter)
          {
            const uint64_t this_hash = std::get<0>(code_pair);
            table->add_hash(this_hash);
          }
        }
      }
    } else
    {
      if (sequence.size() >= qgram_length)
      {
        HashValueIterator qgiter(qgram_length, sequence.data(),
                                 sequence.size());
        for (auto const &&code_pair : qgiter)
        {
          const uint64_t this_hash = std::get<0>(code_pair);
          table->add_hash(this_hash);
        }
      }
    }
  }
}

template <bool split_at_wildcard,
          class HashValueIterator,
          class TableClass>
static TableClass ntcard_enumerate_seq(const std::string &inputfilename,
                                       size_t qgram_length,
                                       size_t s_value,
                                       size_t r_value)
{
  static constexpr const int buf_size = 1 << 14;
  TableClass table(s_value, r_value);
  if (gttl_likely_fasta_format(inputfilename))
  {
    GttlFpType in_fp = gttl_fp_type_open(inputfilename.c_str(), "rb");
    if (in_fp == nullptr)
    {
      throw std::string(": cannot open file");
      /* check_err.py checked */
    }
    GttlSeqIterator<buf_size> gttl_si(in_fp);
    try
    {
      ntcard_enumerate_inner<split_at_wildcard,
                             GttlSeqIterator<buf_size>,
                             HashValueIterator,
                             TableClass>
                            (&gttl_si, &table, qgram_length);
    }
    catch (std::string &msg)
    {
      gttl_fp_type_close(in_fp);
      throw msg;
    }
    gttl_fp_type_close(in_fp);
  } else
  {
    GttlLineIterator<buf_size> line_iterator(inputfilename.c_str());
    using FastQIterator = GttlFastQIterator<GttlLineIterator<buf_size>>;
    FastQIterator fastq_it(line_iterator);
    ntcard_enumerate_inner<split_at_wildcard,
                           FastQIterator,
                           HashValueIterator,
                           TableClass>
                          (&fastq_it, &table, qgram_length);
  }
  return table;
}

template <bool split_at_wildcard,
          class HashValueIterator,
          class TableClass>
static TableClass ntcard_enumerate_thd(const std::string &inputfilename,
                                       size_t qgram_length,
                                       size_t s_value,
                                       size_t r_value,
                                       size_t num_threads)
{
  assert(num_threads > 1);
  const bool fasta_format = gttl_likely_fasta_format(inputfilename);
  SequencesSplit sequence_parts(num_threads, inputfilename, fasta_format);
  TableClass first_table(s_value,r_value);
  std::vector<TableClass *> other_tables;
  for (size_t thd_num = 1; thd_num < sequence_parts.size(); thd_num++)
  {
    other_tables.push_back(new TableClass(s_value, r_value));
  }
  std::vector<std::thread> threads;
  for (size_t thd_num = 0; thd_num < sequence_parts.size(); thd_num++)
  {
    threads.push_back(std::thread([&first_table, &other_tables,
                                   &sequence_parts,thd_num,qgram_length]
    {
      GttlSeqIterator<0> gttl_si(sequence_parts[thd_num]);
      ntcard_enumerate_inner<split_at_wildcard,
                             GttlSeqIterator<0>,
                             HashValueIterator,
                             TableClass>
                            (&gttl_si,
                             thd_num == 0 ? &first_table
                                          : other_tables[thd_num-1],
                             qgram_length);
    }));
  }
  for (auto &th : threads)
  {
    th.join();
  }
  for (size_t thd_num = 1; thd_num < sequence_parts.size(); thd_num++)
  {
    first_table.merge(*other_tables[thd_num-1]);
    delete other_tables[thd_num-1];
  }
  return first_table;
}

template <bool split_at_wildcard,
          class HashValueIterator,
          class TableClass>
static void process_encoded_sequence_part(
  const DNAEncodingMultiLength<uint64_t,split_at_wildcard,false>
    &dna_encoding_multi_length,
  size_t qgram_length,
  size_t part_idx,
  TableClass *table)
{
  const auto end_it = dna_encoding_multi_length.end();
  for(auto it = dna_encoding_multi_length.begin(part_idx); it != end_it; ++it)
  {
    const uint64_t *sub_unit_ptr;
    size_t sequence_length;
    std::tie(sub_unit_ptr, sequence_length) = *it;
    if (sequence_length >= qgram_length)
    {
      std::basic_string<uint8_t> ds = dna_sequence_decode(sub_unit_ptr,
                                                          sequence_length);
      assert(sequence_length == ds.size());
      HashValueIterator qgiter(qgram_length, ds.data(), sequence_length);
      for (auto const &&code_pair : qgiter)
      {
        const uint64_t this_hash = std::get<0>(code_pair);
        table->add_hash(this_hash);
      }
    }
  }
}

template <bool split_at_wildcard,
          class HashValueIterator,
          class TableClass>
static TableClass ntcard_enumerate_thd_gz(const DNAEncodingMultiLength
                                            <uint64_t,split_at_wildcard,false>
                                            &dna_encoding_multi_length,
                                          size_t qgram_length,
                                          size_t s_value,
                                          size_t r_value)
{
  TableClass first_table(s_value,r_value);
  if (dna_encoding_multi_length.num_parts_get() == 1)
  {
    process_encoded_sequence_part<split_at_wildcard,
                                  HashValueIterator,
                                  TableClass>
                                 (dna_encoding_multi_length,
                                  qgram_length,
                                  0,&first_table);
    return first_table;
  }
  std::vector<TableClass *> other_tables;
  for (size_t thd_num = 1; thd_num < dna_encoding_multi_length.num_parts_get();
       thd_num++)
  {
    other_tables.push_back(new TableClass(s_value, r_value));
  }
  std::vector<std::thread> threads;
  for (size_t thd_num = 0;
       thd_num < dna_encoding_multi_length.num_parts_get();
       thd_num++)
  {
    threads.push_back(std::thread([&first_table, &other_tables,
                                   &dna_encoding_multi_length,thd_num,
                                   qgram_length]
    {
      process_encoded_sequence_part<split_at_wildcard,
                                    HashValueIterator,
                                    TableClass>
                                   (dna_encoding_multi_length,
                                    qgram_length,
                                    thd_num,
                                    thd_num == 0 ? &first_table
                                                 : other_tables[thd_num-1]);
    }));
  }
  for (auto &th : threads)
  {
    th.join();
  }
  for (size_t thd_num = 1; thd_num < dna_encoding_multi_length.num_parts_get();
       thd_num++)
  {
    first_table.merge(*other_tables[thd_num-1]);
    delete other_tables[thd_num-1];
  }
  return first_table;
}

template <bool split_at_wildcard,
          class HashValueIterator,
          class HashValueIteratorNoTransform,
          class TableClass>
static TableClass ntcard_enumerate(const std::string &inputfilename,
                                   size_t qgram_length,
                                   size_t s_value,
                                   size_t r_value,
                                   size_t num_threads)
{
  if (num_threads == 1)
  {
    RunTimeClass rt_enumerate{};
    auto table = ntcard_enumerate_seq<split_at_wildcard,
                                      HashValueIterator,
                                      TableClass>
                                     (inputfilename,
                                      qgram_length,
                                      s_value,
                                      r_value);
    rt_enumerate.show("ntcard.enumerate");
    return table;
  }
  if (gttl_has_suffix(inputfilename,std::string(".gz")))
  {
    RunTimeClass rt_encode{};
    DNAEncodingMultiLength<uint64_t,split_at_wildcard,false>
      dna_encoding_multi_length(inputfilename);
    dna_encoding_multi_length.prepare_view(num_threads);
    rt_encode.show("ntcard.encode");
    RunTimeClass rt_enumerate{};
    auto table = ntcard_enumerate_thd_gz<split_at_wildcard,
                                         HashValueIteratorNoTransform,
                                         TableClass>
                                        (dna_encoding_multi_length,
                                         qgram_length,
                                         s_value,
                                         r_value);
    rt_enumerate.show("ntcard.enumerate");
    return table;
  }
  RunTimeClass rt_enumerate{};
  auto table = ntcard_enumerate_thd<split_at_wildcard,
                                    HashValueIterator,
                                    TableClass>
                                   (inputfilename,
                                    qgram_length,
                                    s_value,
                                    r_value,
                                    num_threads);
  rt_enumerate.show("ntcard.enumerate");
  return table;
}
#endif
