/* created by ./score_class_choices.py --alignment_output DO NOT EDIT */
#ifndef ALIGNMENT_OUTPUT_FUNCTION_HPP
#define ALIGNMENT_OUTPUT_FUNCTION_HPP
#include <cstddef>
#include "alignment/blosum62.hpp"
#include "alignment/unit_score_aa.hpp"
#include "alignment/unit_score_nuc.hpp"
#include "alignment/unit_score_nuc_2_2.hpp"
#include "alignment/unit_score_nuc_lower.hpp"
#include "alignment/unit_score_nuc_upper.hpp"
#include "alignment/score_class_base.hpp"
#include "alignment/score_matrix_name.hpp"
#include "sequences/alignment_output.hpp"
#include "sequences/gttl_substring.hpp"

static auto alignment_output_function_get(
                             const char *score_matrix_id,
                             const ScoreMatrixName &score_matrix_name,
                             bool dna_alphabet)
{
  if (!dna_alphabet)
  {
    if (score_matrix_name.is(Score_matrix_undefined)
        || score_matrix_name.is(Score_matrix_blosum62))
    {
      return alignment_output<GttlSubstring<uint8_t>,
                              GttlSubstring<uint8_t>,
                              uint8_t,
                              encoded_matching_characters<Blosum62>,
                              to_char_map<Blosum62>>;
    } else
    {
      if (score_matrix_name.is(Score_matrix_unit_score_aa))
      {
        return alignment_output<GttlSubstring<uint8_t>,
                                GttlSubstring<uint8_t>,
                                uint8_t,
                                encoded_matching_characters<Unit_score_aa>,
                                to_char_map<Unit_score_aa>>;
      } else
      {
        ScoreMatrixName score_matrix_name_instance{};
        StrFormat msg(": score matrix %s is not possible for protein "
                      "sequences; the following choices are available: %s",
                      score_matrix_id,
                      score_matrix_name_instance.string_values_joined(", ")
                                                   .c_str());
        throw std::runtime_error(msg.str());
      }
    }
  } else
  {
    if (score_matrix_name.is(Score_matrix_undefined)
        || score_matrix_name.is(Score_matrix_unit_score_nuc))
    {
      return alignment_output<GttlSubstring<uint8_t>,
                              GttlSubstring<uint8_t>,
                              uint8_t,
                              encoded_matching_characters<Unit_score_nuc>,
                              to_char_map<Unit_score_nuc>>;
    } else
    {
      if (score_matrix_name.is(Score_matrix_undefined)
          || score_matrix_name.is(Score_matrix_unit_score_nuc_2_2))
      {
        return alignment_output<GttlSubstring<uint8_t>,
                                GttlSubstring<uint8_t>,
                                uint8_t,
                                encoded_matching_characters<Unit_score_nuc_2_2>,
                                to_char_map<Unit_score_nuc_2_2>>;
      } else
      {
        if (score_matrix_name.is(Score_matrix_unit_score_nuc_lower))
        {
          return alignment_output<GttlSubstring<uint8_t>,
                                  GttlSubstring<uint8_t>,
                                  uint8_t,
                                  encoded_matching_characters
                                    <Unit_score_nuc_lower>,
                                  to_char_map<Unit_score_nuc_lower>>;
        } else
        {
          if (score_matrix_name.is(Score_matrix_unit_score_nuc_upper))
          {
            return alignment_output<GttlSubstring<uint8_t>,
                                    GttlSubstring<uint8_t>,
                                    uint8_t,
                                    encoded_matching_characters
                                      <Unit_score_nuc_upper>,
                                    to_char_map<Unit_score_nuc_upper>>;
          } else
          {
            ScoreMatrixName score_matrix_name_instance{};
            StrFormat msg(": score matrix %s is not possible for DNA "
                          "sequences; the following choices are "
                          "available: %s",
                          score_matrix_id,
                          score_matrix_name_instance.string_values_joined(", ")
                                                       .c_str());
            throw std::runtime_error(msg.str());
          }
        }
      }
    }
  }
}
#endif /* ALIGNMENT_OUTPUT_FUNCTION_HPP */
