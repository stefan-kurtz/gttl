#ifndef UNIT_SCORE_NUC_HPP
#define UNIT_SCORE_NUC_HPP
#include <cstdint>
struct Unit_score_nuc
{
  static constexpr const char characters[] = "ACGTN";
  static constexpr size_t num_of_chars = size_t(5);
  static constexpr const char character_spec[]
    = "Aa|Cc|Gg|TtUu|NSYWRKVBDHMnsywrkvbdhm";
  static constexpr const int8_t score_matrix[5][5] = {
            /* A   C   G   T   N */
      /* A */ {2, -1, -1, -1, -1},
      /* C */ {-1, 2, -1, -1, -1},
      /* G */ {-1, -1, 2, -1, -1},
      /* T */ {-1, -1, -1, 2, -1},
      /* N */ {-1, -1, -1, -1, -1}};
  static constexpr const int8_t smallest_score = -1;
};
#endif

