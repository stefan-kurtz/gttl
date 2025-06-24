#ifndef SCORE_CLASS_BASE_HPP
#define SCORE_CLASS_BASE_HPP
#include <cassert>
#include <cstddef>
#include <cstdint>
#include "sequences/gttl_multiseq.hpp"
#include "sequences/literate_multiseq.hpp"

template<int alphasize>
static inline int8_t **scorematrix2D_get(const int8_t
                                           matrix[alphasize][alphasize])
{
  int8_t **scorematrix2D = new int8_t * [alphasize];
  scorematrix2D[0] = new int8_t [alphasize * alphasize];
  for (size_t a = 1; a < alphasize; a++)
  {
    scorematrix2D[a] = scorematrix2D[a-1] + alphasize;
  }
  for (size_t a = 0; a < alphasize; a++)
  {
   for (size_t b = 0; b < alphasize; b++)
    {
      scorematrix2D[a][b] = matrix[a][b];
    }
  }
  return scorematrix2D;
}

template<class ScoreClass>
static inline bool encoded_matching_characters(uint8_t a,uint8_t b)
{
  return a < ScoreClass::num_of_chars && a == b;
}

template<class ScoreClass>
static inline char to_char_map(uint8_t cc)
{
  assert(cc <= ScoreClass::num_of_chars);
  return ScoreClass::characters[cc];
}

template<class ScoreClass>
static inline size_t literate_multiseqs(GttlMultiseq *db_multiseq,
                                        GttlMultiseq *query_multiseq)
{
  using LiterateMultiseqScoreClass
    = LiterateMultiseq<ScoreClass::character_spec,ScoreClass::num_of_chars>;
  LiterateMultiseqScoreClass literate_db_multiseq(db_multiseq);
  literate_db_multiseq.perform_sequence_encoding();
  if (query_multiseq != db_multiseq)
  {
    LiterateMultiseqScoreClass literate_query_multiseq(query_multiseq);
    literate_query_multiseq.perform_sequence_encoding();
  }
  return ScoreClass::num_of_chars;
}
#endif
