#ifndef SW_OUTPUT_RESULT_HPP
#define SW_OUTPUT_RESULT_HPP
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <string>
#include <tuple>
#include <string_view>
#include "alignment/score_matrix_name.hpp"
#include <exception>
#include "threading/threads_output_files.hpp"
#include "sequences/gttl_multiseq.hpp"
#include "sequences/eoplist.hpp"
#include "sequences/alignment_output.hpp"
#include "sequences/gttl_substring.hpp"
#include "alignment/blast_stat.hpp"
#include "alignment/loc_align_coords.hpp"
#include "alignment/affine_dband.hpp"
#include "alignment_output_function.hpp"
#include "alignment_display.hpp"

#include "sequence_decode_function.hpp"

static inline bool adpr_simple_matching_chars(uint8_t a,uint8_t b)
{
  return a == b;
}

template<typename ScoreType>
class SWOutputResultShared
{
  using SWAlignmentSequenceInfo
    = AlignmentSequenceInfo<GttlSubstring<uint8_t>,
                            GttlSubstring<uint8_t>>;
  using AlignmentOutputFunction = void (*)(const SWAlignmentSequenceInfo &,
                                           const Eoplist &,
                                           size_t,
                                           size_t,
                                           size_t,
                                           bool,
                                           bool,
                                           size_t,
                                           bool,
                                           size_t,
                                           size_t,
                                           FILE *);
  using SequenceDecodeFunction
    = std::string (*)(const GttlSubstring<uint8_t> &);
  size_t alphasize;
  const int8_t *const *scorematrix2D;
  int8_t smallest_score;
  int8_t gap_open_penalty;
  int8_t gap_extension_penalty;
  uint32_t min_bit_or_raw_score;
  const GttlMultiseq *db_multiseq,
                     *query_multiseq;
  const BlastStatistics *blast_statistics;
  const AlignmentDisplay &alignment_display;
  bool dna_alphabet,
       header_display,
       opt_memory,
       stop_after_first;
  AlignmentOutputFunction alignment_output_function;
  SequenceDecodeFunction sequence_decode_function;
  void write_seqid_or_number(FILE *fpout,size_t i,size_t j) const
  {
    if (header_display)
    {
      assert(db_multiseq != nullptr && query_multiseq != nullptr);
      size_t sh_offset;
      size_t sh_len;
      std::tie(sh_offset,sh_len) = db_multiseq->short_header_get(i);
      const std::string_view db_header = db_multiseq->header_get(i);
      std::fwrite(db_header.data() + sh_offset,sizeof(char),sh_len,fpout);
      std::fputc('\t',fpout);
      std::tie(sh_offset,sh_len) = query_multiseq->short_header_get(j);
      const std::string_view query_header = query_multiseq->header_get(j);
      std::fwrite(query_header.data() + sh_offset,sizeof(char),sh_len,
                  fpout);
      std::fputc('\t',fpout);
    } else
    {
      fprintf(fpout,"%zu\t%zu\t",i,j);
    }
  }

  public:
  SWOutputResultShared(size_t _alphasize,
                       const char *score_matrix_id,
                       ScoreMatrixName score_matrix_name,
                       const int8_t *const *_scorematrix2D,
                       int8_t _smallest_score,
                       int8_t _gap_open_penalty,
                       int8_t _gap_extension_penalty,
                       uint32_t _min_bit_or_raw_score,
                       const GttlMultiseq *_db_multiseq,
                       const GttlMultiseq *_query_multiseq,
                       const BlastStatistics *_blast_statistics,
                       const AlignmentDisplay &_alignment_display,
                       bool _dna_alphabet,
                       bool _header_display,
                       bool _opt_memory,
                       bool _stop_after_first)
    : alphasize(_alphasize)
    , scorematrix2D(_scorematrix2D)
    , smallest_score(_smallest_score)
    , gap_open_penalty(_gap_open_penalty)
    , gap_extension_penalty(_gap_extension_penalty)
    , min_bit_or_raw_score(_min_bit_or_raw_score)
    , db_multiseq(_db_multiseq)
    , query_multiseq(_query_multiseq)
    , blast_statistics(_blast_statistics)
    , alignment_display(_alignment_display)
    , dna_alphabet(_dna_alphabet)
    , header_display(_header_display)
    , opt_memory(_opt_memory)
    , stop_after_first(_stop_after_first)
    , alignment_output_function(nullptr)
    , sequence_decode_function(nullptr)
  {
    try
    {
      alignment_output_function
        = alignment_output_function_get(score_matrix_id,
                                        score_matrix_name,
                                        dna_alphabet);
      sequence_decode_function
        = sequence_decode_function_get(score_matrix_id,
                                       score_matrix_name,
                                       dna_alphabet);
    }
    catch (const std::exception &err)
    {
      assert(false); /* This should not happen as reason for exception must
                        have been detected in scoring_info_and_seq_trans */
    }
  }
  bool process(FILE *fpout, const LocalAlignmentCoordinates &best_coords,
               size_t i, size_t j) const
  {
    const uint32_t score
      = blast_statistics != nullptr
          ? blast_statistics->raw_score2bit_score(best_coords.raw_score)
          : best_coords.raw_score;
    if (score >= min_bit_or_raw_score)
    {
      write_seqid_or_number(fpout,i,j);
      best_coords.show(fpout,dna_alphabet,query_multiseq,j);
    } else
    {
      return false;
    }
    if (blast_statistics != nullptr)
    {
      fprintf(fpout,"\t%u",score);
    }
    if (alignment_display.s_coverage())
    {
      const size_t db_len = db_multiseq->sequence_length_get(i);
      fprintf(fpout,"\t%.2f",100.0 * best_coords.usubstringlength/db_len);
    }
    if (alignment_display.q_coverage())
    {
      const size_t query_len = query_multiseq->sequence_length_get(j);
      fprintf(fpout,"\t%.2f",100.0 * best_coords.vsubstringlength/query_len);
    }
    if (alignment_display.need_alignment() or
        alignment_display.s_substring() or
        alignment_display.q_substring())
    {
      const uint8_t *db_seq_ptr = db_multiseq->encoded_sequence_ptr_get(i);
      const size_t db_len = db_multiseq->sequence_length_get(i);
      const GttlSubstring<uint8_t> usubstring(db_seq_ptr,
                                              best_coords.ustart,
                                              best_coords.usubstringlength);
      const uint8_t *query_seq_ptr
        = query_multiseq->encoded_sequence_ptr_get(j);
      const size_t query_len = query_multiseq->sequence_length_get(j);
      const GttlSubstring<uint8_t> vsubstring(best_coords.forward_strand,
                                              query_seq_ptr,
                                              best_coords.vstart,
                                              best_coords.vsubstringlength,
                                              query_len);
      const bool distinguish_mismatch_match = true;
      Eoplist *eoplist = alignment_display.need_traceback()
                           ? new Eoplist(distinguish_mismatch_match)
                           : nullptr;
      if (alignment_display.need_alignment())
      {
        const bool no_score_run = false;
        assert(best_coords.raw_score > 0);
        GttlAffineDPbanded<ScoreType> adpr(opt_memory,
                                           alignment_display.need_traceback(),
                                           db_len,
                                           query_len);
        const size_t dp_score
          = adpr.template alignment_get<uint8_t,adpr_simple_matching_chars>
                                       (alphasize,
                                        eoplist,
                                        gap_open_penalty,
                                        gap_extension_penalty,
                                        scorematrix2D,
                                        smallest_score,
                                        usubstring,
                                        vsubstring,
                                        no_score_run,
                                        best_coords.raw_score);
        if (dp_score != static_cast<size_t>(best_coords.raw_score))
        {
          fprintf(stderr,"file %s, line %d: dp_score = %zu != %u "
                         "expected_score\n",
                 __FILE__,__LINE__,dp_score,best_coords.raw_score);
          exit(EXIT_FAILURE);
        }
      }
      if (alignment_display.only_verify_score())
      {
        fputc('\n',fpout);
      } else
      {
        if (alignment_display.identity())
        {
          assert(eoplist != nullptr);
          fprintf(fpout,"\t%.2f",100.0 - eoplist->error_percentage_get());
        }
        if (alignment_display.cigar())
        {
          assert(eoplist != nullptr);
          std::string cigar_string
            = eoplist->cigar_string_get(distinguish_mismatch_match);
          fprintf(fpout,"\t%s",cigar_string.c_str());
        }
        if (alignment_display.s_substring())
        {
          std::string usubstring_dec = sequence_decode_function(usubstring);
          fprintf(fpout,"\t%s",usubstring_dec.c_str());
        }
        if (alignment_display.q_substring())
        {
          std::string vsubstring_dec = sequence_decode_function(vsubstring);
          fprintf(fpout,"\t%s",vsubstring_dec.c_str());
        }
        fputc('\n',fpout);
        if (alignment_display.subject_query_alignment())
        {
          SWAlignmentSequenceInfo
            alignment_sequence_info(&usubstring,
                                    best_coords.ustart,
                                    &vsubstring,
                                    best_coords.vstart,
                                    db_multiseq->sequences_maximum_length_get(),
                                    query_multiseq
                                      ->sequences_maximum_length_get());
          static constexpr const size_t top_seqlength = 0;
          static constexpr const size_t low_reference = 0;
          static constexpr const size_t one_off = 0;
          static constexpr const bool subject_first = true;
          static constexpr const bool distinguish_mismatch_match = true;
          alignment_output_function(alignment_sequence_info,
                                    *eoplist,
                                    top_seqlength,
                                    low_reference,
                                    one_off,
                                    subject_first,
                                    distinguish_mismatch_match,
                                    alignment_display.width(),
                                    best_coords.forward_strand,
                                    query_len,
                                    best_coords.vsubstringlength,
                                    fpout);
        }
      }
#ifndef NDEBUG
      if (eoplist != nullptr)
      {
        const auto alignment_score
          = eoplist->evaluate_score<ScoreType,GttlSubstring<uint8_t>,
                                    GttlSubstring<uint8_t>>
                                   (usubstring,
                                    vsubstring,
                                    gap_open_penalty,
                                    gap_extension_penalty,
                                    scorematrix2D);
        if (alignment_score != static_cast<ScoreType>(best_coords.raw_score))
        {
          fprintf(stderr,"file %s, line %d: score = %d != %u "
                         "expected_score\n",
                 __FILE__,__LINE__,alignment_score,best_coords.raw_score);
          exit(EXIT_FAILURE);
        }
      }
#endif
      delete eoplist;
    } else
    {
      fprintf(fpout,"\n");
    }
    return stop_after_first;
  }
};

class SWOutputResultGetThreadRelated
{
  static constexpr const char *progname_prefix = "sw_all_against_all";
  ThreadsOutputFiles *threads_output_files;
  public:
  SWOutputResultGetThreadRelated(size_t num_threads,
                                 const char *threads_out_prefix)
    : threads_output_files(num_threads == 1 ? nullptr
                                            : new ThreadsOutputFiles(
                                                     progname_prefix,
                                                     threads_out_prefix,
                                                     num_threads))
  {}
  FILE *operator [](size_t t) const
  {
    return threads_output_files == nullptr
             ? stdout
             : threads_output_files->filepointer(t);
  }
  ~SWOutputResultGetThreadRelated(void)
  {
    delete threads_output_files;
  }
};
#endif
