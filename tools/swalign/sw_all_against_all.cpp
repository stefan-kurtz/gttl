/*
  Copyright (c) 2022 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2022 Center for Bioinformatics, University of Hamburg

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

#include <cstdint>
#include <exception>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <unordered_map>
#include <tuple>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <utility>
#include <algorithm>
#include <iostream>
#include <cinttypes>
#include <vector>
#include "utilities/matrix_partition.hpp"
#include "utilities/runtime_class.hpp"
#include "utilities/all_vs_all2.hpp"
#include "utilities/gttl_line_iterator.hpp"
#include "utilities/split_string.hpp"
#include "utilities/string_of_digits.hpp"
#include "threading/thread_pool_var.hpp"
#include "sequences/gttl_multiseq.hpp"
#include "alignment/blast_stat.hpp"
#include "sw_option_parser.hpp"
#include "sw_input_data.hpp"
#include "sw_comparator.hpp"
#include "sw_output_result.hpp"
#include "sw_store_best.hpp"

class Restrict2Pairs
{
  using HeaderID2idx = std::unordered_map<std::string,uint32_t>;
  uint32_t convert_header(HeaderID2idx &header_id2idx,const std::string &s)
    const
  {
    if (string_of_digits(s))
    {
      uint32_t number;
      [[maybe_unused]] const int ret =
                                   std::sscanf(s.c_str(), "%" PRIu32, &number);
      assert(ret == 1);
      return number;
    }
    return header_id2idx[s];
  }
  struct HashPair
  {
    size_t operator () (std::pair<uint32_t, uint32_t> const &pair) const
    {
      const size_t combined = (static_cast<size_t>(std::get<0>(pair)) << 32) |
                              static_cast<size_t>(std::get<1>(pair));
      return std::hash<size_t>()(combined);
    }
  };
  using SetOfPairs = std::unordered_set<std::pair<size_t,size_t>,HashPair>;
  SetOfPairs pairs;
  HeaderID2idx extract_header_ids(const GttlMultiseq *multiseq) const
  {
    HeaderID2idx header_id2idx{};
    assert(multiseq->sequences_number_get() > 0 and
           multiseq->sequences_number_get() - 1 <= UINT32_MAX);
    for (uint32_t seqnum = 0;
         seqnum < static_cast<uint32_t>(multiseq->sequences_number_get());
         seqnum++)
    {
      size_t sh_offset;
      size_t sh_len;
      std::tie(sh_offset,sh_len) = multiseq->short_header_get(seqnum);
      const std::string_view header = multiseq->header_get(seqnum);
      const std::string key(header.data() + sh_offset, sh_len);
      header_id2idx[key] = seqnum;
    }
    return header_id2idx;
  }
  SetOfPairs file2pairs_from_indexes(HeaderID2idx &db_header_id2idx,
                                     HeaderID2idx &query_header_id2idx,
                                     const char *inputfile) const
  {
    SetOfPairs local_pairs{};
    constexpr const int buf_size = 1 << 14;
    GttlLineIterator<buf_size> gttl_li(inputfile);
    std::string buffer;
    while (gttl_li.next(&buffer))
    {
      if (buffer.size() > 0 and buffer[0] != '#')
      {
        std::vector<std::string> vec = gttl_split_string(buffer, '\t');
        assert(vec.size() >= 2);
        const uint32_t i = convert_header(db_header_id2idx,vec[0]);
        const uint32_t j = convert_header(query_header_id2idx,vec[1]);
        local_pairs.insert(std::make_pair(i,j));
      }
      buffer.clear();
    }
    if (local_pairs.size() == 0)
    {
      throw std::invalid_argument(
        "file specified with option -r cannot be empty");
    }
    std::cout << "# local_pairs\t" << local_pairs.size() << '\n';
    return local_pairs;
  }
  SetOfPairs filename2pairs_from_multiseqs(const GttlMultiseq *db_multiseq,
                                           const GttlMultiseq *query_multiseq,
                                           const char *inputfile) const
  {
    if (inputfile != nullptr)
    {
      HeaderID2idx db_header_id2idx = extract_header_ids(db_multiseq);
      if (db_multiseq == query_multiseq)
      {
        return file2pairs_from_indexes(db_header_id2idx,db_header_id2idx,
                                       inputfile);
      }
      HeaderID2idx query_header_id2idx
        = extract_header_ids(query_multiseq);
      return file2pairs_from_indexes(db_header_id2idx,query_header_id2idx,
                                     inputfile);
    }
    SetOfPairs empty_local_pairs{};
    return empty_local_pairs;
  }
  public:
  Restrict2Pairs(const GttlMultiseq *db_multiseq,
                 const GttlMultiseq *query_multiseq,const char *inputfile)
    : pairs(filename2pairs_from_multiseqs(db_multiseq,query_multiseq,inputfile))
  {}
  bool check(size_t i, size_t j) const noexcept
  {
    return pairs.size() == 0 or pairs.count(std::make_pair(i,j)) > 0;
  }
};

template<class SWProcessResultShared,class SWProcessResultThreadRelated,
         class SWProcessResultGetThreadRelated>
static void sw_all_against_all(const SWOptions &options,
                               const int8_t *const *scorematrix2D,
                               int8_t smallest_score,
                               size_t alphasize,
                               size_t seqnum_divisor,
                               const GttlMultiseq *db_multiseq,
                               const GttlMultiseq *query_multiseq,
                               const SWProcessResultShared
                                     &sw_process_result_shared,
                               const SWProcessResultGetThreadRelated
                                     &sw_process_result_get_thread_related)
{
  using ThisSWcomparator = SWcomparator<Restrict2Pairs,SWProcessResultShared,
                                        SWProcessResultThreadRelated>;
  std::vector<ThisSWcomparator *> comparator_vector{};
  comparator_vector.reserve(options.num_threads);
  const size_t max_seq_len
    = std::max(db_multiseq->sequences_maximum_length_get(),
               query_multiseq->sequences_maximum_length_get());
  const Restrict2Pairs restrict2pairs(
                               db_multiseq,
                               query_multiseq,
                               options.restrict_to_pairs_file);
  for (size_t t = 0; t < options.num_threads; t++)
  {
    comparator_vector.push_back(new ThisSWcomparator
                                    (sw_process_result_get_thread_related[t],
                                     alphasize,
                                     not options.no_reverse_strand,
                                     scorematrix2D,
                                     smallest_score,
                                     options.gap_open_penalty,
                                     options.gap_extension_penalty,
                                     max_seq_len,
                                     options.vectorized_alignment == 1,
                                     restrict2pairs,
                                     sw_process_result_shared));
  }
  const size_t db_sequences_number = db_multiseq->sequences_number_get();
  const size_t query_sequences_number = query_multiseq->sequences_number_get();
  const size_t cutlen = std::max(size_t(1),
                                 std::max(db_sequences_number,
                                          query_sequences_number)
                                   /seqnum_divisor);
  const MatrixPartition mp = db_multiseq == query_multiseq ?
                               MatrixPartition(cutlen, db_sequences_number) :
                               MatrixPartition(cutlen,
                                               db_sequences_number,
                                               query_sequences_number);
  GttlThreadPoolVar(options.num_threads,
                    mp.size(),
                    all_against_all_compare_pairs<ThisSWcomparator,
                                                  GttlMultiseq>,
                    *db_multiseq,
                    *query_multiseq,
                    db_multiseq == query_multiseq,
                    mp,
                    comparator_vector);
  for (size_t t = 0; t < options.num_threads; t++)
  {
    delete comparator_vector[t];
  }
}

class BestFromThread
{
  const SWResultVector *result;
  public:
  SWResultVector::ConstIterator it;
  BestFromThread(const SWResultVector *_result)
    : result(_result)
    , it(_result->begin())
  {
    assert(it != _result->end());
  }
  bool operator < (const BestFromThread &other) const
  {
    return *it > *other.it;
  }
  bool at_end(void) const
  {
    return it + 1 == result->end();
  }
};

template<typename ScoreType>
static void multiway_merge_results(
   size_t num_threads,
   size_t best,
   const SWOutputResultShared<ScoreType> &sw_output_result_shared,
   const SWStoreBestResultsGetThreadRelated
     &sw_store_best_results_get_thread_related)
{
  std::vector<BestFromThread> all_best{};
  for (size_t t = 0; t < num_threads; t++)
  {
    SWResultVector *result_t = sw_store_best_results_get_thread_related[t];
    result_t->sort();
    all_best.push_back(BestFromThread(result_t));
  }
  for (size_t b = 0; b < best and not all_best.empty(); b++)
  {
    const std::vector<BestFromThread>::iterator next_best = std::min_element(
                                 all_best.begin(), all_best.end());
    const StoredLocalAlignmentCoordinates &slac = *(next_best->it);
    sw_output_result_shared.process(stdout,slac.coords,
                                    slac.u_seqnum,slac.v_seqnum);
    if (next_best->at_end())
    {
      all_best.erase(next_best);
    } else
    {
      ++next_best->it;
    }
  }
}

int main(int argc,char *argv[])
{
  const char *program_name = argv[0];
  bool haserr = false;
  SWOptions options{true};
  if (argc == 1)
  {
    options.usage(program_name);
    return EXIT_SUCCESS;
  }
  try
  {
    options.parse(argc, argv);
  }
  catch (const std::exception &err)
  {
    std::cerr << program_name << ": " << err.what() << '\n';
    haserr = true;
  }
  RunTimeClass timer{};

  GttlMultiseq *db_multiseq = nullptr;
  GttlMultiseq *query_multiseq = nullptr;
  bool dna_alphabet = false;
  int8_t **scorematrix2D = nullptr;
  int8_t smallest_score = INT8_MAX;
  size_t alphasize = 0;
  BlastStatistics *blast_statistics = nullptr;

  if (!haserr)
  {
    try
    {
      std::tie(db_multiseq,query_multiseq,dna_alphabet,
               scorematrix2D,smallest_score,alphasize,blast_statistics)
        = sw_input_data(options);
    }
    catch(const std::exception &err)
    {
      std::cerr << program_name << ": file \"" << options.dbfile << "\""
                << err.what() << '\n';
      haserr = true;
    }
  }
  if (not dna_alphabet)
  {
    options.no_reverse_strand = true;
  }
  if (!haserr)
  {
    options.show(stdout,program_name);
    options.show_fields(stdout,
                        dna_alphabet,
                        blast_statistics != nullptr);

    using ScoreType = int32_t;
    const SWOutputResultShared<ScoreType> sw_output_result_shared(
                                 alphasize,
                                 options.score_matrix_id,
                                 options.score_matrix_name,
                                 static_cast<const int8_t *const *>(
                                                              scorematrix2D),
                                 smallest_score,
                                 options.gap_open_penalty,
                                 options.gap_extension_penalty,
                                 options.min_bit_score,
                                 db_multiseq,
                                 query_multiseq,
                                 blast_statistics,
                                 options.alignment_display,
                                 dna_alphabet,
                                 options.header_display,
                                 options.opt_memory,
                                 options.stop_after_first);
    const size_t seqnum_divisor = options.stop_after_first ? size_t(1)
                                                           : size_t(10);
    try
    {
      if (options.best == 0)
      {
        const SWOutputResultGetThreadRelated
          sw_output_result_get_thread_related(
                                     options.num_threads,
                                     options.threads_out_prefix);

        sw_all_against_all<SWOutputResultShared<ScoreType>,FILE,
                           SWOutputResultGetThreadRelated>
                          (options,
                           static_cast<const int8_t *const *>(scorematrix2D),
                           smallest_score,
                           alphasize,
                           seqnum_divisor,
                           db_multiseq,
                           query_multiseq,
                           sw_output_result_shared,
                           sw_output_result_get_thread_related);
      } else
      {
        const SWStoreBestResultsShared sw_store_best_results_shared{};
        const SWStoreBestResultsGetThreadRelated
                                     sw_store_best_results_get_thread_related(
                                                        options.num_threads,
                                                        options.best);
        sw_all_against_all<SWStoreBestResultsShared,SWResultVector,
                           SWStoreBestResultsGetThreadRelated>
                          (options,
                           static_cast<const int8_t *const *>(scorematrix2D),
                           smallest_score,
                           alphasize,
                           seqnum_divisor,
                           db_multiseq,
                           query_multiseq,
                           sw_store_best_results_shared,
                           sw_store_best_results_get_thread_related);
        if (options.num_threads == 1)
        {
          SWResultVector *result_0
            = sw_store_best_results_get_thread_related[0];
          result_0->sort();
          for (auto &&slac : *result_0)
          {
            sw_output_result_shared.process(stdout,slac.coords,
                                            slac.u_seqnum,slac.v_seqnum);
          }
        } else
        {
          multiway_merge_results<ScoreType>
                                (options.num_threads,
                                 options.best,
                                 sw_output_result_shared,
                                 sw_store_best_results_get_thread_related);
        }
      }
    }
    catch (const std::exception &err)
    {
      std::cerr << program_name << ": " << err.what() << '\n';
      haserr = true;
    }
  }
  if (scorematrix2D != nullptr)
  {
    delete[] scorematrix2D[0];
    delete[] scorematrix2D;
  }
  if (db_multiseq != query_multiseq)
  {
    delete query_multiseq;
  }
  delete db_multiseq;
  delete blast_statistics;
  //ksw_show_score_matrix_accesses();
  if (!haserr)
  {
    timer.show("sw_all_against_all");
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
