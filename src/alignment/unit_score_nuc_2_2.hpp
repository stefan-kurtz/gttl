#ifndef UNIT_SCORE_NUC_2_2_HPP
#define UNIT_SCORE_NUC_2_2_HPP
#include <cstdint>
/* score matrix corresponding to affine cost function
   0,4,6.2 */
struct Unit_score_nuc_2_2
{
  static constexpr const char characters[] = "ACGTN";
  static constexpr size_t num_of_chars = size_t(5);
  static constexpr const char character_spec[]
    = "Aa|Cc|Gg|TtUu|NSYWRKVBDHMnsywrkvbdhm";
  static constexpr const int8_t score_matrix[5][5] = {
            /* A   C   G   T   N */
      /* A */ {2, -2, -2, -2, -2},
      /* C */ {-2, 2, -2, -2, -2},
      /* G */ {-2, -2, 2, -2, -2},
      /* T */ {-2, -2, -2, 2, -2},
      /* N */ {-2, -2, -2, -2, -2}};
  static constexpr const int8_t smallest_score = -2;
};
#endif
