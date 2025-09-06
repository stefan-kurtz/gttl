#ifndef NTCARD_HPP
#define NTCARD_HPP

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <ios>
#include <string>
#include <string_view>
#include <thread>
#include <tuple>
#include <vector>
#include <format>
#include "sequences/gttl_fasta_generator.hpp"
#include "sequences/gttl_fastq_generator.hpp"
#include "utilities/gttl_file_open.hpp"
#include "utilities/has_fasta_or_fastq_extension.hpp"
#include "sequences/split.hpp"
#include "sequences/dna_seq_encoder.hpp"
#include "sequences/dna_seq_decoder.hpp"
#include "utilities/runtime_class.hpp"

template <bool split_at_wildcard,
          class SeqGenerator,
          class HashValueIterator,
          class TableClass,
          bool is_aminoacid>
static void ntcard_enumerate_inner(SeqGenerator* gttl_si,
                                   TableClass* table,
                                   size_t qgram_length)
{
  size_t sequences_number = 0;
  for (auto &&si : *gttl_si)
  {
    auto sequence = si->sequence_get();
    if constexpr (split_at_wildcard)
    {
      const NtCardRanger<is_aminoacid> nuc_ranger(
                                   sequence.data(), sequence.size());
      for (auto const &&range : nuc_ranger)
      {
        const size_t this_length = std::get<1>(range);
        if (this_length >= qgram_length)
        {
          const char *const substring = sequence.data() + std::get<0>(range);
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
    sequences_number++;
  }
  table->sequences_number_set(sequences_number);
}

template <bool split_at_wildcard,
          class HashValueIterator,
          class TableClass,
          bool is_aminoacid>
static TableClass ntcard_enumerate_seq(const std::string &inputfilename,
                                       size_t qgram_length,
                                       size_t s_value,
                                       size_t r_value)
{
  static constexpr const int buf_size = 1 << 14;
  TableClass table(s_value, r_value);
  if (gttl_likely_fasta_format(inputfilename))
  {
    const GttlFpType in_fp = gttl_fp_type_open(inputfilename.c_str(), "rb");
    if (in_fp == nullptr)
    {
      throw std::ios_base::failure(": cannot open file");
      /* check_err.py checked */
    }
    GttlFastAGenerator<buf_size> gttl_si(in_fp);
    ntcard_enumerate_inner<split_at_wildcard,
                           GttlFastAGenerator<buf_size>,
                           HashValueIterator,
                           TableClass,
                           is_aminoacid>
                          (&gttl_si, &table, qgram_length);
  } else
  {
    GttlFastQGenerator<buf_size> fastq_it(inputfilename.c_str());
    ntcard_enumerate_inner<split_at_wildcard,
                           GttlFastQGenerator<buf_size>,
                           HashValueIterator,
                           TableClass,
                           is_aminoacid>
                          (&fastq_it, &table, qgram_length);
  }
  return table;
}

template <bool split_at_wildcard,
          class HashValueIterator,
          class TableClass,
          bool is_aminoacid>
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
  constexpr const size_t buf_size = size_t{1} << size_t{14};
  if (gttl_likely_fasta_format(inputfilename))
  {
    for (size_t thd_num = 0; thd_num < sequence_parts.size(); thd_num++)
    {
      threads.push_back(std::thread([&first_table, &other_tables,
                                     &sequence_parts,thd_num,qgram_length]
      {
        GttlFastAGenerator<buf_size> gttl_si(sequence_parts[thd_num]);
        ntcard_enumerate_inner<split_at_wildcard,
                               GttlFastAGenerator<buf_size>,
                               HashValueIterator,
                               TableClass,
                               is_aminoacid>
                              (&gttl_si,
                               thd_num == 0 ? &first_table
                                            : other_tables[thd_num-1],
                               qgram_length);
      }));
    }
  } else
  {
    for (size_t thd_num = 0; thd_num < sequence_parts.size(); thd_num++)
    {
      threads.push_back(std::thread([&first_table, &other_tables,
                                     &sequence_parts,thd_num,qgram_length]
      {
        const std::string_view &this_view =  sequence_parts[thd_num];
        GttlFastQGenerator<buf_size> fastq_it(this_view.data(),
                                              this_view.size());
        ntcard_enumerate_inner<split_at_wildcard,
                               GttlFastQGenerator<buf_size>,
                               HashValueIterator,
                               TableClass,
                               is_aminoacid>
                              (&fastq_it,
                               thd_num == 0 ? &first_table
                                            : other_tables[thd_num-1],
                               qgram_length);
      }));
    }
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
          class TableClass,
          bool is_aminoacid>
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

  threads.reserve(dna_encoding_multi_length.num_parts_get());
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
          class TableClass,
          bool is_aminoacid>
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
                                      TableClass,
                                      is_aminoacid>
                                     (inputfilename,
                                      qgram_length,
                                      s_value,
                                      r_value);
    const std::string msg = std::format("ntcard.enumerate, {}, {}, 1 thread",
                                        inputfilename,
                                        is_aminoacid ? "protein" : "DNA");
    rt_enumerate.show(msg.c_str());
    return table;
  }
  if (inputfilename.ends_with(".gz"))
  {
    RunTimeClass rt_encode{};
    DNAEncodingMultiLength<uint64_t,split_at_wildcard,false>
      dna_encoding_multi_length(inputfilename);
    dna_encoding_multi_length.prepare_view(num_threads);
    rt_encode.show("ntcard.encode");
    RunTimeClass rt_enumerate{};
    auto table = ntcard_enumerate_thd_gz<split_at_wildcard,
                                         HashValueIteratorNoTransform,
                                         TableClass,
                                         is_aminoacid>
                                        (dna_encoding_multi_length,
                                         qgram_length,
                                         s_value,
                                         r_value);
    table.sequences_number_set(dna_encoding_multi_length
                                  .total_number_of_sequences_get());
    const std::string msg = std::format("ntcard.enumerate {}, {}, {} threads",
                                        inputfilename,
                                        is_aminoacid ? "protein" : "DNA",
                                        num_threads);
    rt_enumerate.show(msg.c_str());
    return table;
  }
  RunTimeClass rt_enumerate{};
  auto table = ntcard_enumerate_thd<split_at_wildcard,
                                    HashValueIterator,
                                    TableClass,
                                    is_aminoacid>
                                   (inputfilename,
                                    qgram_length,
                                    s_value,
                                    r_value,
                                    num_threads);
  const std::string msg = std::format("ntcard.enumerate, {}, {}, {} threads",
                                      inputfilename,
                                      is_aminoacid ? "protein" : "DNA",
                                      num_threads);
  rt_enumerate.show(msg.c_str());
  return table;
}
#endif
