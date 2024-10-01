#ifndef ALIGNMENT_OUTPUT_HPP
#define ALIGNMENT_OUTPUT_HPP
#include <cmath>
#include <cstdio>
#include <cstddef>
#include <cmath>
#include <cassert>
#include <iostream>
#include "utilities/unused.hpp"
#include "sequences/eoplist.hpp"

template<class SeqClass_u,class SeqClass_v>
struct AlignmentSequenceInfo
{
  const SeqClass_u *useq;
  size_t ustart;
  const SeqClass_v *vseq;
  size_t vstart;
  size_t max_length;
  AlignmentSequenceInfo (void)
    : useq(nullptr)
    , ustart(0)
    , vseq(nullptr)
    , vstart(0)
    , max_length(0)
  {}

  AlignmentSequenceInfo(const SeqClass_u *_useq,
                        size_t _ustart,
                        const SeqClass_v *_vseq,
                        size_t _vstart,
                        size_t umax_length,
                        size_t vmax_length)
    : useq(_useq)
    , ustart(_ustart)
    , vseq(_vseq)
    , vstart(_vstart)
    , max_length(std::max(umax_length,vmax_length))
  {}
  void set(const SeqClass_u *_useq,
           size_t _ustart,
           const SeqClass_v *_vseq,
           size_t _vstart,
           size_t ulen,
           size_t vlen)
  {
    useq = _useq;
    vseq = _vseq;
    ustart = _ustart;
    vstart = _vstart;
    max_length = std::max(ulen,vlen);
  }
  int width_of_numbers_get(void) const
  {
    assert(max_length > 0);
    const size_t maximum_position = max_length-1;
    if (maximum_position < 10)
    {
      return 1;
    }
    return 1 + std::log10(static_cast<double>(maximum_position));
  }
};

static void alignment_output_single_line(const char *tag,
                                         int width_of_numbers,
                                         size_t width_alignment,
                                         const char *buf,
                                         size_t start_pos,
                                         size_t end_pos,
                                         FILE *fp)
{
  std::fprintf(fp,"%s  %-*zu  ",tag,width_of_numbers,start_pos);
  std::fwrite(buf,sizeof *buf,width_alignment,fp);
  std::fprintf(fp,"  %zu\n",end_pos);
}

static void alignment_output_middle_line(int width_of_numbers,
                                         size_t width_alignment,
                                         const char *midbuf,
                                         FILE *fp)
{
  static constexpr const int width_of_Sbjct_Query = 5;
  /* 5 is the length of the strings Sbjct and Query */
  std::fprintf(fp,"%*s",width_of_numbers + width_of_Sbjct_Query + 4,"");
  std::fwrite(midbuf,sizeof *midbuf,width_alignment,fp);
  std::fputc('\n',fp);
}

static void alignment_output_write_lines(size_t one_off,
                                         bool subject_first,
                                         int width_of_numbers,
                                         size_t width_alignment,
                                         const char *subject_buf,
                                         size_t subject_seqlength,
                                         size_t subject_start_pos,
                                         size_t subject_end_pos,
                                         const char *midbuf,
                                         const char *query_buf,
                                         size_t query_start_pos,
                                         size_t query_end_pos,
                                         FILE *fp)
{
  assert(width_of_numbers > 0);
  if (subject_first)
  {
    alignment_output_single_line("Sbjct",
                                 width_of_numbers,
                                 width_alignment,
                                 subject_buf,
                                 one_off + subject_start_pos,
                                 one_off + subject_end_pos,
                                 fp);
    alignment_output_middle_line(width_of_numbers,width_alignment,midbuf,fp);
    alignment_output_single_line("Query",
                                 width_of_numbers,
                                 width_alignment,
                                 query_buf,
                                 one_off + query_start_pos,
                                 one_off + query_end_pos,
                                 fp);
  } else
  {
    alignment_output_single_line("Query",
                                 width_of_numbers,
                                 width_alignment,
                                 query_buf,
                                 query_start_pos + one_off,
                                 query_end_pos + one_off,
                                 fp);
    alignment_output_middle_line(width_of_numbers,width_alignment,midbuf,fp);
    if (subject_seqlength == 0)
    {
      alignment_output_single_line("Sbjct",
                                   width_of_numbers,
                                   width_alignment,
                                   subject_buf,
                                   subject_start_pos + one_off,
                                   subject_end_pos + one_off,fp);
    } else
    {
      assert(subject_seqlength > subject_start_pos &&
                subject_seqlength >= subject_end_pos);
      alignment_output_single_line("Sbjct",
                                   width_of_numbers,
                                   width_alignment,
                                   subject_buf,
                                   subject_seqlength - 1
                                     - subject_start_pos + one_off,
                                   one_off + (subject_seqlength >
                                              subject_end_pos
                                                 ? subject_seqlength - 1
                                                     - subject_end_pos
                                                 : 0),
                                   fp);
    }
  }
  fputc('\n',fp);
}

static inline size_t alignment_output_show_advance(size_t one_off,
                                            bool subject_first,
                                            int width_of_numbers,
                                            size_t pos,
                                            size_t width_alignment,
                                            const char *topbuf,
                                            size_t top_seqlength,
                                            size_t top_start_pos,
                                            size_t top_end_pos,
                                            const char *midbuf,
                                            const char *lowbuf,
                                            size_t low_start_pos,
                                            size_t low_end_pos,
                                            FILE *fp)
{
  assert(width_alignment > 0);
  if (pos + 1 < width_alignment)
  {
    return pos + 1;
  }
  assert(pos == width_alignment - 1);
  alignment_output_write_lines(one_off,
                               subject_first,
                               width_of_numbers,
                               width_alignment,
                               topbuf,
                               top_seqlength,
                               top_start_pos,
                               top_end_pos,
                               midbuf,
                               lowbuf,
                               low_start_pos,
                               low_end_pos,
                               fp);
  return 0;
}

template<class SeqClass_u,
         class SeqClass_v,
         typename CharType,
         bool (&match_method)(CharType,CharType),
         char (&to_char)(CharType)>
static void alignment_output(const AlignmentSequenceInfo <SeqClass_u,SeqClass_v>
                               &asi,
                             const Eoplist &eoplist,
                             size_t top_seqlength,
                             size_t low_reference,
                             size_t one_off,
                             bool subject_first,
                             GTTL_UNUSED bool distinguish_mismatch_match,
                             size_t width_alignment,
                             bool forward_strand,
                             size_t low_seq_length,
                             size_t low_match_length,
                             FILE *fp)
{
  static constexpr const char alignment_output_matchsymbol = '|';
  static constexpr const char alignment_output_mismatchsymbol = ' ';
  static constexpr const char alignment_output_gapsymbol = '-';
  const int width_of_numbers = asi.width_of_numbers_get();
  const size_t low_start_base
    = forward_strand ? (low_reference == 0 ? asi.vstart
                                           : (low_reference - asi.vstart))
                     : (low_seq_length - asi.vstart - low_match_length);
  size_t pos = 0,
         idx_u = 0,
         idx_v = 0,
         top_start_pos = asi.ustart,
         low_start_pos = low_start_base;

  char *topbuf = new char [3 * width_alignment],
       *midbuf = topbuf + width_alignment,
       *lowbuf = midbuf + width_alignment;

  const SeqClass_u &useq = *(asi.useq);
  const SeqClass_v &vseq = *(asi.vseq);
  for (auto &&co : eoplist)
  {
    switch (co.edit_operation)
    {
      case MatchOp:
      case MismatchOp:
        for (size_t j = 0; j < co.iteration && idx_u < useq.size() &&
                                               idx_v < vseq.size(); j++)
        {
          char cc_a, cc_b;
          bool is_match;

          is_match = match_method(useq[idx_u],vseq[idx_v]);
          cc_a = to_char(useq[idx_u]);
          cc_b = to_char(vseq[idx_v]);
          assert(pos < width_alignment);
          topbuf[pos] = cc_a;
          lowbuf[pos] = cc_b;
          if (distinguish_mismatch_match)
          {
            if (co.edit_operation == MatchOp)
            {
              assert(is_match);
            } else
            {
              assert(!is_match);
            }
          }
          midbuf[pos] = is_match ? alignment_output_matchsymbol
                                 : alignment_output_mismatchsymbol;
          pos = alignment_output_show_advance(one_off,
                                              subject_first,
                                              width_of_numbers,
                                              pos,
                                              width_alignment,
                                              topbuf,
                                              top_seqlength,
                                              top_start_pos,
                                              asi.ustart + idx_u,
                                              midbuf,
                                              lowbuf,
                                              low_start_pos,
                                              low_start_base + idx_v,
                                              fp);
          if (pos == 0)
          {
            top_start_pos = asi.ustart + idx_u + 1;
            low_start_pos = low_start_base + idx_v + 1;
          }
          idx_u++;
          idx_v++;
        }
        break;
      case DeletionOp:
        for (size_t j = 0; j < co.iteration && idx_u < useq.size(); j++)
        {
          const char cc_a = to_char(useq[idx_u]);
          assert(pos < width_alignment);
          topbuf[pos] = cc_a;
          midbuf[pos] = alignment_output_mismatchsymbol;
          lowbuf[pos] = alignment_output_gapsymbol;
          pos = alignment_output_show_advance(one_off,
                                              subject_first,
                                              width_of_numbers,
                                              pos,
                                              width_alignment,
                                              topbuf,
                                              top_seqlength,
                                              top_start_pos,
                                              asi.ustart + idx_u,
                                              midbuf,
                                              lowbuf,
                                              low_start_pos,
                                              low_start_base + idx_v - 1,
                                              fp);
          if (pos == 0)
          {
            top_start_pos = asi.ustart + idx_u + 1;
            low_start_pos = low_start_base + idx_v;
          }
          idx_u++;
        }
        break;
      case InsertionOp:
        for (size_t j = 0; j < co.iteration && idx_v < vseq.size(); j++)
        {
          const char cc_b = to_char(vseq[idx_v]);
          assert(pos < width_alignment);
          topbuf[pos] = alignment_output_gapsymbol;
          midbuf[pos] = alignment_output_mismatchsymbol;
          lowbuf[pos] = cc_b;
          pos = alignment_output_show_advance(one_off,
                                              subject_first,
                                              width_of_numbers,
                                              pos,
                                              width_alignment,
                                              topbuf,
                                              top_seqlength,
                                              top_start_pos,
                                              asi.ustart + idx_u - 1,
                                              midbuf,
                                              lowbuf,
                                              low_start_pos,
                                              low_start_base + idx_v,
                                              fp);
          if (pos == 0)
          {
            top_start_pos = asi.ustart + idx_u;
            low_start_pos = low_start_base + idx_v + 1;
          }
          idx_v++;
        }
        break;
      default:
        assert(false);
    }
  }
  if (pos > 0)
  {
    alignment_output_write_lines(one_off,
                                 subject_first,
                                 width_of_numbers,
                                 pos,
                                 topbuf,
                                 top_seqlength,
                                 top_start_pos,
                                 asi.ustart + std::min(idx_u,useq.size() - 1),
                                 midbuf,
                                 lowbuf,
                                 low_start_pos,
                                 low_start_base + std::min(idx_v,
                                                           vseq.size() - 1),
                                 fp);
  }
  delete[] topbuf;
}
#endif
