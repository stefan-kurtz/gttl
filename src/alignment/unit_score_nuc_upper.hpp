/* generated by ./scorematrix.py --constexpr unit_score_nuc_upper
 * unit_score_nuc_upper.txt DO NOT EDIT */
#ifndef UNIT_SCORE_NUC_UPPER_HPP
#define UNIT_SCORE_NUC_UPPER_HPP
#include <cstdint>
#include <cstddef>
struct Unit_score_nuc_upper
{
  static constexpr const char characters[] = "ACGTN";
  static constexpr size_t num_of_chars = size_t(5);
  static constexpr const char character_spec[] = "A|C|G|T|N";
  static constexpr const int8_t score_matrix[5][5] = {
      /* A  C  G  T  N */ /* A */ {2, -1, -1, -1, -1},
      /* C */ {-1, 2, -1, -1, -1}, /* G */ {-1, -1, 2, -1, -1},
      /* T */ {-1, -1, -1, 2, -1}, /* N */ {-1, -1, -1, -1, -1}};
  static constexpr const int8_t smallest_score = -1;
};
#endif
