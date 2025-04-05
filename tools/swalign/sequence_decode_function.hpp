/* created by ./score_class_choices.py --sequence_decode DO NOT EDIT */
#ifndef SEQUENCE_DECODE_FUNCTION_HPP
#define SEQUENCE_DECODE_FUNCTION_HPP
#include <cstddef>
#include "alignment/blosum62.hpp"
#include "alignment/unit_score_aa.hpp"
#include "alignment/unit_score_nuc.hpp"
#include "alignment/unit_score_nuc_2_2.hpp"
#include "alignment/unit_score_nuc_lower.hpp"
#include "alignment/unit_score_nuc_upper.hpp"
#include "alignment/score_class_base.hpp"
#include "alignment/score_matrix_name.hpp"
#include "sequences/gttl_substring.hpp"

template <class SeqClass, typename CharType, char (&to_char)(CharType)>
static std::string sequence_decode(const SeqClass &seq)
{
  std::string s{};
  for (size_t idx = 0; idx < seq.size(); idx++)
  {
    s += to_char(seq[idx]);
  }
  return s;
}

static auto sequence_decode_function_get(
    const char *score_matrix_id, const ScoreMatrixName &score_matrix_name,
    bool dna_alphabet)
{
  if (!dna_alphabet)
  {
    if (score_matrix_name.is(Score_matrix_undefined) ||
        score_matrix_name.is(Score_matrix_blosum62))
    {
      return sequence_decode<GttlSubstring<uint8_t>, uint8_t,
                             to_char_map<Blosum62>>;
    } else
    {
      if (score_matrix_name.is(Score_matrix_unit_score_aa))
      {
        return sequence_decode<GttlSubstring<uint8_t>, uint8_t,
                               to_char_map<Unit_score_aa>>;
      } else
      {
        ScoreMatrixName score_matrix_name_instance{};
        StrFormat msg(
            ": score matrix %s is not possible for protein "
            "sequences; the following choices are available: %s",
            score_matrix_id,
            score_matrix_name_instance.string_values_joined(", ").c_str());
        throw msg.str();
      }
    }
  } else
  {
    if (score_matrix_name.is(Score_matrix_undefined) ||
        score_matrix_name.is(Score_matrix_unit_score_nuc))
    {
      return sequence_decode<GttlSubstring<uint8_t>, uint8_t,
                             to_char_map<Unit_score_nuc>>;
    } else
    {
      if (score_matrix_name.is(Score_matrix_undefined) ||
          score_matrix_name.is(Score_matrix_unit_score_nuc_2_2))
      {
        return sequence_decode<GttlSubstring<uint8_t>, uint8_t,
                               to_char_map<Unit_score_nuc_2_2>>;
      } else
      {
        if (score_matrix_name.is(Score_matrix_unit_score_nuc_lower))
        {
          return sequence_decode<GttlSubstring<uint8_t>, uint8_t,
                                 to_char_map<Unit_score_nuc_lower>>;
        } else
        {
          if (score_matrix_name.is(Score_matrix_unit_score_nuc_upper))
          {
            return sequence_decode<GttlSubstring<uint8_t>, uint8_t,
                                   to_char_map<Unit_score_nuc_upper>>;
          } else
          {
            ScoreMatrixName score_matrix_name_instance{};
            StrFormat msg(
                ": score matrix %s is not possible for DNA "
                "sequences; the following choices are "
                "available: %s",
                score_matrix_id,
                score_matrix_name_instance.string_values_joined(", ").c_str());
            throw msg.str();
          }
        }
      }
    }
  }
}
#endif /* SEQUENCE_DECODE_FUNCTION_HPP */
