/* created by ./score_class_choices.py --scoring DO NOT EDIT */
#ifndef SCORING_INFO_AND_SEQ_TRANS_HPP
#define SCORING_INFO_AND_SEQ_TRANS_HPP
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <tuple>
#include "alignment/blosum62.hpp"
#include "alignment/unit_score_aa.hpp"
#include "alignment/unit_score_nuc.hpp"
#include "alignment/unit_score_nuc_2_2.hpp"
#include "alignment/unit_score_nuc_lower.hpp"
#include "alignment/unit_score_nuc_upper.hpp"
#include "alignment/score_class_base.hpp"
#include "alignment/score_matrix_name.hpp"
#include "sequences/gttl_multiseq.hpp"
#include "utilities/str_format.hpp"

static inline std::tuple<int8_t **, int8_t, size_t> scoring_info_and_seq_trans(
    const char *score_matrix_id, const ScoreMatrixName &score_matrix_name,
    bool dna_alphabet, GttlMultiseq *db_multiseq, GttlMultiseq *query_multiseq)
{
  if (!dna_alphabet)
  {
    if (score_matrix_name.is(Score_matrix_undefined) ||
        score_matrix_name.is(Score_matrix_blosum62))
    {
      return {scorematrix2D_get<Blosum62::num_of_chars>(Blosum62::score_matrix),
              Blosum62::smallest_score,
              literate_multiseqs<Blosum62>(db_multiseq, query_multiseq)};
    } else
    {
      if (score_matrix_name.is(Score_matrix_unit_score_aa))
      {
        return {scorematrix2D_get<Unit_score_aa::num_of_chars>(
                    Unit_score_aa::score_matrix),
                Unit_score_aa::smallest_score,
                literate_multiseqs<Unit_score_aa>(db_multiseq, query_multiseq)};
      } else
      {
        const ScoreMatrixName score_matrix_name_instance{};
        const StrFormat msg(": score matrix %s is not possible for protein "
                            "sequences; the following choices are available:"
                            "%s",
                            score_matrix_id,
                            score_matrix_name_instance
                              .string_values_joined(", ")
                              .c_str());
        throw std::runtime_error(msg.str());
      }
    }
  } else
  {
    if (score_matrix_name.is(Score_matrix_undefined) ||
        score_matrix_name.is(Score_matrix_unit_score_nuc))
    {
      return {scorematrix2D_get<Unit_score_nuc::num_of_chars>(
                  Unit_score_nuc::score_matrix),
              Unit_score_nuc::smallest_score,
              literate_multiseqs<Unit_score_nuc>(db_multiseq, query_multiseq)};
    } else
    {
      if (score_matrix_name.is(Score_matrix_undefined) ||
          score_matrix_name.is(Score_matrix_unit_score_nuc_2_2))
      {
        return {scorematrix2D_get<Unit_score_nuc_2_2::num_of_chars>(
                    Unit_score_nuc_2_2::score_matrix),
                Unit_score_nuc_2_2::smallest_score,
                literate_multiseqs<Unit_score_nuc_2_2>(db_multiseq,
                                                       query_multiseq)};
      } else
      {
        if (score_matrix_name.is(Score_matrix_unit_score_nuc_lower))
        {
          return {scorematrix2D_get<Unit_score_nuc_lower::num_of_chars>(
                      Unit_score_nuc_lower::score_matrix),
                  Unit_score_nuc_lower::smallest_score,
                  literate_multiseqs<Unit_score_nuc_lower>(db_multiseq,
                                                           query_multiseq)};
        } else
        {
          if (score_matrix_name.is(Score_matrix_unit_score_nuc_upper))
          {
            return {scorematrix2D_get<Unit_score_nuc_upper::num_of_chars>(
                        Unit_score_nuc_upper::score_matrix),
                    Unit_score_nuc_upper::smallest_score,
                    literate_multiseqs<Unit_score_nuc_upper>(db_multiseq,
                                                             query_multiseq)};
          } else
          {
            const ScoreMatrixName score_matrix_name_instance{};
            const StrFormat msg(": score matrix %s is not possible for DNA "
                                "sequences; the following choices are "
                                "available: %s",
                                score_matrix_id,
                                score_matrix_name_instance
                                  .string_values_joined(", ")
                                  .c_str());
            throw std::runtime_error(msg.str());
          }
        }
      }
    }
  }
  return std::tuple<int8_t **, int8_t, size_t>(NULL, INT8_MAX, 0);
}
#endif /* SCORING_INFO_AND_SEQ_TRANS_HPP */
