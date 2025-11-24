#ifndef AFFINE_DBAND_HPP
#define AFFINE_DBAND_HPP
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cstddef>
#include <cassert>
#include <algorithm>
#include <utility>
#include "sequences/eoplist.hpp"
#include "sequences/gttl_substring.hpp"

template<typename ScoreType>
class GttlAffineDPbanded
{
  using AffineAlignEditOp = enum : uint8_t
  {
    AFFINE_DBAND_UNDEF = 0, /* unknown, keep it, as it' the default value 0 */
    AFFINE_DBAND_REPLACEMENT, /* 1 = 01, shift 0 = 2 * (1-1) */
    AFFINE_DBAND_DELETION,    /* 2 = 10, shift 2 = 2 * (2-1) */
    AFFINE_DBAND_INSERTION    /* 3 = 11, shift 4 = 2 * (3-1) */
  };

  using AffineAlignScoreTriple = struct {
    ScoreType Rvalue,
              Dvalue,
              Ivalue;
  };
  /* <AffineAlignTraceBits> objects describe the backtracing edges
     relating to the last edit operation R,D,I. We store all three values
     between 0 and 3 as bitpairs in one uint8_t-value, using the shift
     values as described above. */

  class AffineAlignTraceBits
  {
    uint8_t trace;
    public:
    template<int edge_type>
    [[nodiscard]] AffineAlignEditOp edge_get(void) const noexcept
    {
      constexpr const int shift = 2 * (edge_type-1);
      return static_cast<AffineAlignEditOp>((trace >> shift) & uint8_t(3));
    }
    template<int edge_type>
    void edge_set(void)
    {
      constexpr const int shift = 2 * (edge_type-1);
      trace = static_cast<uint8_t>(edge_type << shift);
    }
    void edge_set_all(AffineAlignEditOp rmaxedge,
                      AffineAlignEditOp dmaxedge,
                      AffineAlignEditOp imaxedge)
    {
      constexpr const int shift_r = 2 * (AFFINE_DBAND_REPLACEMENT-1);
      constexpr const int shift_d = 2 * (AFFINE_DBAND_DELETION-1);
      constexpr const int shift_i = 2 * (AFFINE_DBAND_INSERTION-1);
      trace = (rmaxedge << shift_r) |
              (dmaxedge << shift_d) |
              (imaxedge << shift_i);
    }
  };

  size_t base_type_size,
         max_ulen,
         max_vlen,
         matrix_cells;
  AffineAlignScoreTriple *columnspace;
  void *matrix_column_reservoir,
       *matrix_space;

  template<typename CharType>
  ScoreType affine_diagonalband_fillDPtab_bits(
                                      [[maybe_unused]] size_t alphasize,
                                      AffineAlignTraceBits **bitmatrix,
                                      int8_t gap_opening, /* >= 0 */
                                      int8_t gap_extension, /* > 0 */
                                      const int8_t * const *scorematrix2D,
                                      ScoreType min_align_score,
                                      const GttlSubstring<CharType>
                                        &usubstring,
                                      const GttlSubstring<CharType>
                                        &vsubstring,
                                      int64_t left_dist,
                                      int64_t right_dist)
  {
    const ScoreType start_penalty = gap_opening + gap_extension;
#ifndef NDEBUG
    const size_t band_width = static_cast<size_t>(right_dist - left_dist + 1);
    const int64_t lendiff = static_cast<int64_t>(vsubstring.size()) -
                            static_cast<int64_t>(usubstring.size());
#endif /* NDEBUG */

    assert(bitmatrix != nullptr &&
           columnspace != nullptr &&
           gap_opening >= 0 &&
           gap_extension > 0 &&
           left_dist <= std::min(int64_t(0), lendiff) &&
           left_dist >= -static_cast<int64_t>(usubstring.size()) &&
           right_dist >= std::max(int64_t(0), lendiff) &&
           right_dist <= static_cast<int64_t>(vsubstring.size()));
    size_t high_row = static_cast<size_t>(-left_dist);

    columnspace->Rvalue = 0;
    columnspace->Dvalue = -gap_opening;
    columnspace->Ivalue = -gap_opening;
    AffineAlignTraceBits *colptr = bitmatrix[0];
    for (size_t i = 1; i <= high_row; i++)
    {
      bitmatrix[0][i].template edge_set<AFFINE_DBAND_DELETION>();
      columnspace[i].Rvalue = min_align_score;
      columnspace[i].Dvalue = columnspace[i-1].Dvalue - gap_extension;
      columnspace[i].Ivalue = min_align_score;
    }
    size_t low_row = 0;
    CharType *const usubstring_cache = new CharType[usubstring.size()];
    for (size_t idx = 0; idx < usubstring.size(); idx++)
    {
      usubstring_cache[idx] = usubstring[idx];
    }
    for (size_t j = 1; j <= vsubstring.size(); j++)
    {
      const CharType cb = vsubstring[j-1];
      assert(cb < alphasize);

      ScoreType first_ivalue = min_align_score;
      const size_t prev_high_row = high_row;

      assert(low_row <= high_row &&
             static_cast<size_t>(high_row - low_row + 1) <= band_width);
      colptr += (high_row - low_row + 1);
      bitmatrix[j] = colptr - low_row;
      /* below diagonal band*/
      assert(right_dist >= 0);
      if (std::cmp_less_equal(j, right_dist))
      {
        assert(low_row <= prev_high_row);
        first_ivalue = columnspace[low_row].Ivalue - gap_extension;
        bitmatrix[j][low_row].template edge_set<AFFINE_DBAND_INSERTION>();
      }
      AffineAlignScoreTriple nw = columnspace[low_row];
      columnspace[low_row].Rvalue = min_align_score;
      columnspace[low_row].Dvalue = min_align_score;
      columnspace[low_row].Ivalue = first_ivalue;

      if (high_row < usubstring.size())
      {
        high_row++;
      }
      const int8_t *const score_row = scorematrix2D[static_cast<size_t>(cb)];
      /* diagonalband */
      for (size_t i = low_row+1; i <= high_row; i++)
      {
        AffineAlignScoreTriple currententry;
        AffineAlignEditOp rmaxedge = AFFINE_DBAND_REPLACEMENT;

        /* compute A_affine(i,j,R) from value in the north west */
        currententry.Rvalue = nw.Rvalue;
        if (currententry.Rvalue < nw.Dvalue)
        {
          currententry.Rvalue = nw.Dvalue;
          rmaxedge = AFFINE_DBAND_DELETION;
        }
        if (currententry.Rvalue < nw.Ivalue)
        {
          currententry.Rvalue = nw.Ivalue;
          rmaxedge = AFFINE_DBAND_INSERTION;
        }
        const size_t ca_idx = static_cast<size_t>(usubstring_cache[i-1]);
        currententry.Rvalue += static_cast<ScoreType>(score_row[ca_idx]);
        /* compute A_affine(i,j,D) from value in the north */
        const ScoreType score_from_R = columnspace[i-1].Rvalue - start_penalty;
        const ScoreType score_from_D = columnspace[i-1].Dvalue - gap_extension;
        currententry.Dvalue = score_from_R >= score_from_D ? score_from_R
                                                           : score_from_D;
        const AffineAlignEditOp dmaxedge = score_from_R >= score_from_D ?
                                     AFFINE_DBAND_REPLACEMENT :
                                     AFFINE_DBAND_DELETION;
        currententry.Ivalue = min_align_score;
        AffineAlignEditOp imaxedge = AFFINE_DBAND_UNDEF;
        /* compute A_affine(i,j,I) from value in the west, if available*/
        if (i <= prev_high_row)
        {
          const ScoreType score_from_R = columnspace[i].Rvalue - start_penalty;
          const ScoreType score_from_I = columnspace[i].Ivalue - gap_extension;
          currententry.Ivalue = score_from_R >= score_from_I ? score_from_R
                                                             : score_from_I;
          imaxedge = score_from_R >= score_from_I ? AFFINE_DBAND_REPLACEMENT
                                                  : AFFINE_DBAND_INSERTION;
        }
        nw = columnspace[i];
        columnspace[i] = currententry;
        bitmatrix[j][i].edge_set_all(rmaxedge,dmaxedge,imaxedge);
      }
      assert(low_row < high_row && right_dist >= 0);
      if (std::cmp_greater(j, right_dist))
      {
        low_row++;
      }
    }
    delete[] usubstring_cache;
    return columnspace[usubstring.size()].Rvalue;
  }

  template<bool keep_columns,typename CharType>
  ScoreType affine_diagonalband_fillDPtab_scores(
                                      [[maybe_unused]] size_t alphasize,
                                      [[maybe_unused]] /* if not keep_columns */
                                      AffineAlignScoreTriple **dpmatrix,
                                      int8_t gap_opening, /* >= 0 */
                                      int8_t gap_extension, /* > 0 */
                                      const int8_t * const *scorematrix2D,
                                      ScoreType min_align_score,
                                      const GttlSubstring<CharType>
                                        &usubstring,
                                      const GttlSubstring<CharType>
                                        &vsubstring,
                                      int64_t left_dist,
                                      int64_t right_dist)
  {
    const ScoreType start_penalty = gap_opening + gap_extension;
    [[maybe_unused]] AffineAlignScoreTriple *colptr;
  #ifndef NDEBUG
    const size_t band_width = static_cast<size_t>(right_dist - left_dist + 1);
    const int64_t lendiff = static_cast<int64_t>(vsubstring.size()) -
                            static_cast<int64_t>(usubstring.size());
  #endif /* NDEBUG */

    assert(columnspace != nullptr &&
           gap_opening >= 0 &&
           gap_extension > 0 &&
           left_dist <= std::min(int64_t(0), lendiff) &&
           left_dist >= -static_cast<int64_t>(usubstring.size()) &&
           right_dist >= std::max(int64_t(0), lendiff) &&
           right_dist <= static_cast<int64_t>(vsubstring.size()));
    size_t high_row = static_cast<size_t>(-left_dist);

    /* first entry */
    columnspace->Rvalue = 0;
    columnspace->Dvalue = -start_penalty;
    columnspace->Ivalue = -start_penalty;
    /* first column */
    for (size_t i = 1; i <= high_row; i++)
    {
      assert(usubstring[i-1] < alphasize);
      columnspace[i].Rvalue = min_align_score;
      columnspace[i].Dvalue = columnspace[i-1].Dvalue - gap_extension;
      columnspace[i].Ivalue = min_align_score;
    }
    size_t low_row = 0;
    if constexpr (keep_columns)
    {
      assert(dpmatrix != nullptr);
      const size_t width = high_row - low_row + 1;
      colptr = dpmatrix[0];
      memcpy(colptr,columnspace,width * sizeof *columnspace);
      colptr += width;
    } else
    {
      colptr = nullptr;
    }
    CharType *const usubstring_cache = new CharType[usubstring.size()];
    for (size_t idx = 0; idx < usubstring.size(); idx++)
    {
      usubstring_cache[idx] = usubstring[idx];
    }
    /* next columns */
    for (size_t j = 1; j <= vsubstring.size(); j++)
    {
      const CharType cb = vsubstring[j-1];
      assert(cb < alphasize);
      ScoreType first_ivalue = min_align_score;
      const size_t prev_high_row = high_row;

      assert(low_row <= high_row &&
             static_cast<size_t>(high_row - low_row + 1) <= band_width);
      if (std::cmp_less_equal(j, right_dist))
      {
        assert(low_row <= prev_high_row);
        first_ivalue = columnspace[low_row].Ivalue - gap_extension;
      }
      AffineAlignScoreTriple nw = columnspace[low_row];
      columnspace[low_row].Rvalue = min_align_score;
      columnspace[low_row].Dvalue = min_align_score;
      columnspace[low_row].Ivalue = first_ivalue;

      const int8_t *const score_row = scorematrix2D[static_cast<size_t>(cb)];
      /* do not branch in the inner loop (except for maximum computation) */
      for (size_t i = low_row+1; i <= prev_high_row; i++)
      {
        const size_t ca_idx = static_cast<size_t>(usubstring_cache[i-1]);
        const ScoreType score_from_R1 = columnspace[i-1].Rvalue - start_penalty;
        const ScoreType score_from_D = columnspace[i-1].Dvalue - gap_extension;
        const ScoreType score_from_R2 = columnspace[i].Rvalue - start_penalty;
        const ScoreType score_from_I = columnspace[i].Ivalue - gap_extension;
        AffineAlignScoreTriple currententry;
        currententry.Dvalue = std::max(score_from_R1,score_from_D);
        currententry.Ivalue = std::max(score_from_R2,score_from_I);
        currententry.Rvalue = std::max(nw.Rvalue,std::max(nw.Dvalue,nw.Ivalue))
                              + static_cast<ScoreType>(score_row[ca_idx]);
        nw = columnspace[i];
        columnspace[i] = currententry;
      }
      if (high_row < usubstring.size())
      {
        AffineAlignScoreTriple currententry;
        const size_t ca_idx = static_cast<size_t>(usubstring_cache[high_row]);

        high_row++;
        const ScoreType score_from_R = columnspace[prev_high_row].Rvalue
                                       - start_penalty;
        const ScoreType score_from_D = columnspace[prev_high_row].Dvalue
                                       - gap_extension;
        currententry.Dvalue = std::max(score_from_R,score_from_D);
        currententry.Ivalue = min_align_score;
        currententry.Rvalue = std::max(nw.Rvalue,std::max(nw.Dvalue,nw.Ivalue))
                              + static_cast<ScoreType>(score_row[ca_idx]);
        columnspace[high_row] = currententry;
      }
      if (std::cmp_greater(j, right_dist))
      {
        low_row++;
      }
      if constexpr (keep_columns)
      {
        const size_t width = high_row - low_row + 1;
        dpmatrix[j] = colptr - low_row;
        memcpy(colptr,columnspace + low_row,width * sizeof *columnspace);
        colptr += width;
      }
      assert(low_row < high_row &&
             ((j <= static_cast<size_t>(right_dist) && low_row == 0) ||
              (j > static_cast<size_t>(right_dist) &&
               low_row == j - static_cast<size_t>(right_dist))) &&
             (high_row == std::min(usubstring.size(),
                                   j + static_cast<size_t>(-left_dist))));
    }
    delete[] usubstring_cache;
    return columnspace[usubstring.size()].Rvalue;
  }

  template<typename CharType,bool (&match_method)(CharType,CharType)>
  void affine_global_alignment_traceback_bits(
                                  Eoplist *eoplist,
                                  const AffineAlignTraceBits * const *bitmatrix,
                                  const GttlSubstring<CharType> &usubstring,
                                  const GttlSubstring<CharType> &vsubstring)
  {
    assert(eoplist != nullptr && bitmatrix != nullptr);
    AffineAlignEditOp edge = AFFINE_DBAND_REPLACEMENT;
    size_t i = usubstring.size();
    size_t j = vsubstring.size();
    while (i > 0 || j > 0)
    {
      const AffineAlignTraceBits trace_bits = bitmatrix[j][i];
      if (edge == AFFINE_DBAND_REPLACEMENT)
      {
        assert(i > 0 && j > 0);
        if (match_method(usubstring[i-1],vsubstring[j-1]))
        {
          eoplist->match_add(1);
        } else
        {
          eoplist->mismatch_add();
        }
        edge = trace_bits.template edge_get<AFFINE_DBAND_REPLACEMENT>();
        i--;
        j--;
      } else
      {
        if (edge == AFFINE_DBAND_DELETION)
        {
          eoplist->deletion_add();
          edge = trace_bits.template edge_get<AFFINE_DBAND_DELETION>();
          assert(i > 0);
          i--;
        } else
        {
          assert(edge == AFFINE_DBAND_INSERTION);
          eoplist->insertion_add();
          edge = trace_bits.template edge_get<AFFINE_DBAND_INSERTION>();
          assert(j > 0);
          j--;
        }
      }
    }
  }

#ifndef NDEBUG
  size_t low_row_on_the_fly(size_t j,int64_t right_dist)
  {
    assert(right_dist >= 0);
    return j <= static_cast<size_t>(right_dist)
             ? 0
             : j - static_cast<size_t>(right_dist);
  }

  size_t high_row_on_the_fly(size_t j,size_t usubstringlen,
                                    int64_t left_dist)
  {
    assert(left_dist <= 0);
    return std::min(usubstringlen,j + static_cast<size_t>(-left_dist));
  }
#endif /* NDEBUG */

  template<typename CharType,bool (&match_method)(CharType,CharType)>
  void affine_global_alignment_traceback_scores(
                                      Eoplist *eoplist,
                                      [[maybe_unused]] /* if not keep_columns */
                                      const AffineAlignScoreTriple *
                                        const *dpmatrix,
                                      int8_t gap_opening, /* > 0 */
                                      int8_t gap_extension, /* > 0 */
                                      const GttlSubstring<CharType>
                                        &usubstring,
                                      const GttlSubstring<CharType>
                                        &vsubstring,
                                      [[maybe_unused]] int64_t left_dist,
                                      [[maybe_unused]] int64_t right_dist)
  {
    const ScoreType start_penalty = gap_opening + gap_extension;

    assert(eoplist != nullptr && dpmatrix != nullptr);
    AffineAlignEditOp edge = AFFINE_DBAND_REPLACEMENT;
    const AffineAlignScoreTriple *previous;
    size_t i = usubstring.size();
    size_t j = vsubstring.size();
    while (i > 0 || j > 0)
    {
      if (edge == AFFINE_DBAND_REPLACEMENT)
      {
        assert(i > 0 && j > 0);
        if (match_method(usubstring[i-1],vsubstring[j-1]))
        {
          eoplist->match_add(1);
        } else
        {
          eoplist->mismatch_add();
        }
        j--;
        i--;
        assert(i >= low_row_on_the_fly(j,right_dist) &&
               i <= high_row_on_the_fly(j,usubstring.size(),left_dist));
        previous = dpmatrix[j] + i;
        const ScoreType maxvalue = std::max(previous->Rvalue,
                                            std::max(previous->Dvalue,
                                                     previous->Ivalue));
        if (maxvalue > previous->Rvalue)
        {
          if (maxvalue == previous->Dvalue)
          {
            edge = AFFINE_DBAND_DELETION;
          } else
          {
            if (maxvalue == previous->Ivalue)
            {
              edge = AFFINE_DBAND_INSERTION;
            }
          }
        }
      } else
      {
        if (edge == AFFINE_DBAND_DELETION)
        {
          eoplist->deletion_add();
          i--;
          assert(i >= low_row_on_the_fly(j,right_dist) &&
                 i <= high_row_on_the_fly(j,usubstring.size(),left_dist));
          previous = dpmatrix[j] + i;
          if (previous->Rvalue - start_penalty >=
              previous->Dvalue - gap_extension)
          {
            edge = AFFINE_DBAND_REPLACEMENT;
          }
        } else
        {
          assert(edge == AFFINE_DBAND_INSERTION);
          assert(j > 0);
          eoplist->insertion_add();
          j--;
          assert(i >= low_row_on_the_fly(j,right_dist) &&
                 i <= high_row_on_the_fly(j,usubstring.size(),left_dist));
          previous = dpmatrix[j] + i;
          if (previous->Rvalue - start_penalty >=
              previous->Ivalue - gap_extension)
          {
            edge = AFFINE_DBAND_REPLACEMENT;
          }
        }
      }
    }
  }

  [[nodiscard]] size_t next_band_width(size_t band_width) const noexcept
  {
    if (band_width < 4)
    {
      return band_width * 2;
    }
    if (band_width < 20)
    {
      return (band_width * 3)/2;
    }
    return (band_width * 5)/4;
  }

  [[nodiscard]] bool affine_opt_memory(void) const noexcept
  {
    return base_type_size == sizeof(AffineAlignTraceBits);
  }

  template<typename T>
  T *col_ptr_get(size_t vsubstringlen)
  {
    if (vsubstringlen > max_vlen)
    {
      max_vlen = std::max(vsubstringlen,
                          static_cast<size_t>(max_vlen * 1.2 + 128));
      matrix_column_reservoir = realloc(matrix_column_reservoir,
                                        (max_vlen+1) * sizeof (void *));
    } else
    {
      assert(matrix_column_reservoir != nullptr);
    }
    return reinterpret_cast<T *>(matrix_column_reservoir);
  }

  void column_space_adjust(size_t usubstringlen)
  {
    if (usubstringlen > max_ulen)
    {
      max_ulen = std::max(usubstringlen,
                          static_cast<size_t>(max_ulen * 1.2 + 128));
      columnspace = static_cast<AffineAlignScoreTriple *>
                               (realloc(columnspace,
                                        (max_ulen+1) * sizeof *columnspace));
    } else
    {
      assert(columnspace != nullptr);
    }
  }

  template<typename T>
  T *matrix_space_get(size_t band_width, size_t vsubstringlen)
  {
    if (band_width * (vsubstringlen+1) >= matrix_cells)
    {
      matrix_cells = std::max(band_width * (vsubstringlen+1),
                              static_cast<size_t>(matrix_cells * 1.2 + 1024));
      matrix_space = realloc(matrix_space, matrix_cells * base_type_size);
    }
    return reinterpret_cast<T *>(matrix_space);
  }

/* calculate alignment within diagonalband specified by left_dist and
   right_dist.  space and running time is O(bandwidth * vsubstringlen) */
  template<typename CharType>
  ScoreType lastcolumnRvalue_get(size_t alphasize,
                                 bool keep_columns,
                                 int8_t gap_opening, /* >= 0 */
                                 int8_t gap_extension, /* > 0 */
                                 const int8_t * const *scorematrix2D,
                                 ScoreType min_align_score,
                                 const GttlSubstring<CharType> &usubstring,
                                 const GttlSubstring<CharType> &vsubstring,
                                 int64_t left_dist,
                                 int64_t right_dist)
  {
    column_space_adjust(usubstring.size());

    if (keep_columns)
    {
      const size_t band_width = static_cast<size_t>(right_dist - left_dist + 1);
      if (affine_opt_memory())
      {
        AffineAlignTraceBits **const bitmatrix = col_ptr_get<
                                     AffineAlignTraceBits *>(vsubstring.size());

        bitmatrix[0] = matrix_space_get<AffineAlignTraceBits>
                                       (band_width,vsubstring.size());
        return affine_diagonalband_fillDPtab_bits<CharType>
                                                 (alphasize,
                                                  bitmatrix,
                                                  gap_opening,
                                                  gap_extension,
                                                  scorematrix2D,
                                                  min_align_score,
                                                  usubstring,
                                                  vsubstring,
                                                  left_dist,
                                                  right_dist);
      } else
      {
        AffineAlignScoreTriple **const dpmatrix = col_ptr_get<
                                     AffineAlignScoreTriple *>(
                                     vsubstring.size());

        dpmatrix[0] = matrix_space_get<AffineAlignScoreTriple>
                                      (band_width,vsubstring.size());
        return affine_diagonalband_fillDPtab_scores<true,CharType>
                                                   (alphasize,
                                                    dpmatrix,
                                                    gap_opening,
                                                    gap_extension,
                                                    scorematrix2D,
                                                    min_align_score,
                                                    usubstring,
                                                    vsubstring,
                                                    left_dist,
                                                    right_dist);
      }
    } else
    {
      return affine_diagonalband_fillDPtab_scores<false,CharType>
                                                 (alphasize,
                                                  nullptr,
                                                  gap_opening,
                                                  gap_extension,
                                                  scorematrix2D,
                                                  min_align_score,
                                                  usubstring,
                                                  vsubstring,
                                                  left_dist,
                                                  right_dist);
    }
  }

  template<typename CharType,bool (&match_method)(CharType,CharType)>
  void traceback(Eoplist *eoplist,
                 int8_t gap_opening, /* >= 0 */
                 int8_t gap_extension, /* > 0 */
                 const GttlSubstring<CharType> &usubstring,
                 const GttlSubstring<CharType> &vsubstring,
                 int64_t left_dist,
                 int64_t right_dist)
  {
    eoplist->reset();
    if (affine_opt_memory())
    {
      const AffineAlignTraceBits *const *const bitmatrix = col_ptr_get<
                                   const AffineAlignTraceBits *const>(
                                   vsubstring.size());
      affine_global_alignment_traceback_bits<CharType,match_method>
                                            (eoplist,
                                             bitmatrix,
                                             usubstring,
                                             vsubstring);
    } else
    {
      const AffineAlignScoreTriple *const *const dpmatrix = col_ptr_get<
                                   const AffineAlignScoreTriple *const>(
                                   vsubstring.size());
      affine_global_alignment_traceback_scores<CharType,match_method>
                                              (eoplist,
                                               dpmatrix,
                                               gap_opening,
                                               gap_extension,
                                               usubstring,
                                               vsubstring,
                                               left_dist,
                                               right_dist);
    }
    eoplist->reverse_end(0);
  }

  public:
  GttlAffineDPbanded (bool opt_memory,
                         bool need_eoplist,
                         size_t _max_ulen,
                         size_t _max_vlen)
    : base_type_size(opt_memory ? sizeof(AffineAlignTraceBits)
                                : sizeof(AffineAlignScoreTriple))
    , max_ulen(_max_ulen)
    , max_vlen(_max_vlen)
    , matrix_cells(0)
    , columnspace(_max_ulen > 0
                    ? static_cast<AffineAlignScoreTriple *>
                      (malloc((max_ulen+1) * sizeof *columnspace))
                    : nullptr)
    , matrix_column_reservoir(need_eoplist && _max_vlen > 0
                                ? malloc((max_vlen+1) * sizeof (void *))
                                : nullptr)
    , matrix_space(nullptr)
  {
  }
  ~GttlAffineDPbanded (void)
  {
    free(matrix_column_reservoir);
    free(matrix_space);
    free(columnspace);
  }
  template<typename CharType,bool (&match_method)(CharType,CharType)>
  size_t alignment_get(size_t alphasize,
                       Eoplist *eoplist,
                       int8_t gap_opening, /* >= 0 */
                       int8_t gap_extension, /* > 0 */
                       const int8_t * const *scorematrix2D,
                       int8_t smallest_score,
                       const GttlSubstring<CharType> &usubstring,
                       const GttlSubstring<CharType> &vsubstring,
                       bool no_score_run,
                       size_t expected_score)
  {
    //std::cout << "usubstring=" << usubstring.to_string() << std::endl;
    //std::cout << "vsubstring=" << vsubstring.to_string() << std::endl;
#ifndef NDEBUG
    const int64_t lendiff = static_cast<int64_t>(vsubstring.size()) -
                            static_cast<int64_t>(usubstring.size());
#endif /* NDEBUG */
    size_t band_width = 1 + static_cast<size_t>
                                       (vsubstring.size() > usubstring.size()
                                          ? vsubstring.size()
                                            - usubstring.size()
                                          : usubstring.size()
                                            - vsubstring.size());
    const ScoreType min_align_score
      = (usubstring.size() + vsubstring.size())
        * static_cast<ScoreType>(smallest_score);
    ScoreType previous_dpscore = min_align_score;

    assert(smallest_score < 0 && min_align_score > INT_MIN/2);
    for (/*Nothing*/; /* Nothing */; /*Nothing*/)
    {
      int64_t left_dist = -static_cast<int64_t>(band_width);
      int64_t right_dist = static_cast<int64_t>(band_width);

      assert(left_dist <= std::min(int64_t(0),lendiff) &&
             right_dist >= std::max(int64_t(0),lendiff));
      left_dist = std::max(left_dist, -static_cast<int64_t>(usubstring.size()));
      right_dist = std::min(right_dist,
                            static_cast<int64_t>(vsubstring.size()));
      const ScoreType dpscore
        = lastcolumnRvalue_get<CharType>(alphasize,
                                         eoplist != nullptr,
                                         gap_opening,
                                         gap_extension,
                                         scorematrix2D,
                                         static_cast<ScoreType>
                                                    (min_align_score),
                                         usubstring,
                                         vsubstring,
                                         left_dist,
                                         right_dist);
      if (expected_score == 0 ||
          dpscore >= static_cast<ScoreType>(expected_score) ||
          (no_score_run && previous_dpscore == dpscore))
      {
        if (eoplist != nullptr)
        {
          traceback<CharType,match_method>(eoplist,
                                           gap_opening, /* >= 0 */
                                           gap_extension, /* > 0 */
                                           usubstring,
                                           vsubstring,
                                           left_dist,
                                           right_dist);
        }
        return static_cast<size_t>(dpscore);
      }
      previous_dpscore = dpscore;
      band_width = next_band_width(band_width);
      assert (left_dist != -static_cast<int64_t>(usubstring.size()) ||
              right_dist != static_cast<int64_t>(vsubstring.size()));
    }
    assert(false);
    return 0;
  }
};
#endif
