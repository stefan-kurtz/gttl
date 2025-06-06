#ifndef SW_INPUT_DATA_HPP
#define SW_INPUT_DATA_HPP
#include <exception>
#include "sequences/gttl_multiseq.hpp"
#include "sequences/multiseq_pair.hpp"
#include "alignment/blast_stat.hpp"
#include "alignment/scoring_info_and_seq_trans.hpp"
#include "sw_option_parser.hpp"

static auto sw_input_data(const SWOptions &options)
{
  GttlMultiseq *db_multiseq = nullptr;
  GttlMultiseq *query_multiseq = nullptr;
  bool dna_alphabet = false;
  std::tie(db_multiseq,query_multiseq,dna_alphabet)
    = create_multiseq_pair(options.dbfile,options.queryfile);

  if (options.header_display or options.restrict_to_pairs_file != nullptr)
  {
    assert(db_multiseq != NULL);
    db_multiseq->short_header_cache_create<'|','|'>();
    if (db_multiseq != query_multiseq)
    {
      assert(query_multiseq != NULL);
      query_multiseq->short_header_cache_create<'|','|'>();
    }
  }

  int8_t **scorematrix2D = NULL;
  int8_t smallest_score = INT8_MAX;
  size_t alphasize = 0;

  std::tie(scorematrix2D,smallest_score,alphasize)
    = scoring_info_and_seq_trans(options.score_matrix_id,
                                 options.score_matrix_name,
                                 dna_alphabet,
                                 db_multiseq,
                                 query_multiseq);

  BlastStatistics *blast_stat = NULL;
  if (!dna_alphabet && (options.score_matrix_name.is(Score_matrix_undefined) ||
                        options.score_matrix_name.is(Score_matrix_blosum62)))
  {
    try
    {
      const bool scaled = false;
      blast_stat = new BlastStatistics(options.gap_open_penalty,
                                       options.gap_extension_penalty,
                                       scaled);
    }
    catch (const std::exception &err)
    {
      if (options.min_bit_score > 0)
      {
        throw err;
      }
    }
  }
  return std::tuple<GttlMultiseq *,GttlMultiseq *,bool,
                    int8_t **,int8_t,size_t,BlastStatistics *>
          (db_multiseq,query_multiseq,dna_alphabet,
           scorematrix2D,smallest_score,alphasize,blast_stat);
}
#endif
