/* created by fastedist_gencode.py --relative --template 10. DO NOT EDIT. */
#ifndef OUTSENSEEDIST_UNR_HPP
#define OUTSENSEEDIST_UNR_HPP
#include <cstddef>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include "sequences/lcs_lcp_len_type.hpp"
#include "sequences/outsenseedist_inplace.hpp"
#ifndef MAX2
#define MAX2(X, Y) ((X) < (Y) ? (Y) : (X))
#endif
template <typename char_type, typename FrontValue,
          LcsLcpLenType suffix_or_prefix_match_len>
static size_t fastedist_unrolled(size_t d_max, const char_type *useq,
                                 size_t ulen, const char_type *vseq,
                                 size_t vlen, size_t useqnum, size_t vseqnum)
{
  const size_t d_bound = vlen - ulen;
  if (d_max == 0)
  {
    fprintf(stderr, "%s: illegal value d_max=0: can only handle d_max > 0\n",
            __func__);
    exit(EXIT_FAILURE);
  }
  assert(ulen < vlen);
  FrontValue front[40];
  front[0] = suffix_or_prefix_match_len(useq, 0, vseq, 0, ulen, vlen, useqnum,
                                        vseqnum);
  /* d = 1: compute front[1...3] */
  /* h = -1 */
  front[37] = front[0] + 1;
  if (front[37] < ulen)
  {
    front[37] += suffix_or_prefix_match_len(
        useq, front[37], vseq, front[37] - 1, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 0 */
  front[38] = front[0] + 1;
  if (front[38] < ulen)
  {
    front[38] += suffix_or_prefix_match_len(useq, front[38], vseq, front[38],
                                            ulen, vlen, useqnum, vseqnum);
  }
  /* h = 1 */
  front[39] = front[0];
  if (front[39] < ulen && front[39] + 1 < vlen)
  {
    front[39] += suffix_or_prefix_match_len(
        useq, front[39], vseq, front[39] + 1, ulen, vlen, useqnum, vseqnum);
  }
  if (d_bound <= 1 && front[38 + d_bound] >= ulen)
  {
    return 1;
  }
  if (d_max == 1)
  {
    return 1 + 1;
  }
  /* d = 2: compute front[4...8] */
  /* h = -2 */
  front[0] = front[37] + 1;
  if (front[0] < ulen)
  {
    front[0] += suffix_or_prefix_match_len(useq, front[0], vseq, front[0] - 2,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -1 */
  front[1] = MAX2(front[37], front[38]) + 1;
  if (front[1] < ulen)
  {
    front[1] += suffix_or_prefix_match_len(useq, front[1], vseq, front[1] - 1,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 0 */
  front[2] = MAX2(front[37], MAX2(front[38], front[39]) + 1);
  if (front[2] < ulen)
  {
    front[2] += suffix_or_prefix_match_len(useq, front[2], vseq, front[2], ulen,
                                           vlen, useqnum, vseqnum);
  }
  /* h = 1 */
  front[3] = MAX2(front[38], front[39] + 1);
  if (front[3] < ulen && front[3] + 1 < vlen)
  {
    front[3] += suffix_or_prefix_match_len(useq, front[3], vseq, front[3] + 1,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 2 */
  front[4] = front[39];
  if (front[4] < ulen && front[4] + 2 < vlen)
  {
    front[4] += suffix_or_prefix_match_len(useq, front[4], vseq, front[4] + 2,
                                           ulen, vlen, useqnum, vseqnum);
  }
  if (d_bound <= 2 && front[2 + d_bound] >= ulen)
  {
    return 2;
  }
  if (d_max == 2)
  {
    return 2 + 1;
  }
  /* d = 3: compute front[9...15] */
  /* h = -3 */
  front[33] = front[0] + 1;
  if (front[33] < ulen)
  {
    front[33] += suffix_or_prefix_match_len(
        useq, front[33], vseq, front[33] - 3, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -2 */
  front[34] = MAX2(front[0], front[1]) + 1;
  if (front[34] < ulen)
  {
    front[34] += suffix_or_prefix_match_len(
        useq, front[34], vseq, front[34] - 2, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -1 */
  front[35] = MAX2(front[0], MAX2(front[1], front[2]) + 1);
  if (front[35] < ulen)
  {
    front[35] += suffix_or_prefix_match_len(
        useq, front[35], vseq, front[35] - 1, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 0 */
  front[36] = MAX2(front[1], MAX2(front[2], front[3]) + 1);
  if (front[36] < ulen)
  {
    front[36] += suffix_or_prefix_match_len(useq, front[36], vseq, front[36],
                                            ulen, vlen, useqnum, vseqnum);
  }
  /* h = 1 */
  front[37] = MAX2(front[2], MAX2(front[3], front[4]) + 1);
  if (front[37] < ulen && front[37] + 1 < vlen)
  {
    front[37] += suffix_or_prefix_match_len(
        useq, front[37], vseq, front[37] + 1, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 2 */
  front[38] = MAX2(front[3], front[4] + 1);
  if (front[38] < ulen && front[38] + 2 < vlen)
  {
    front[38] += suffix_or_prefix_match_len(
        useq, front[38], vseq, front[38] + 2, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 3 */
  front[39] = front[4];
  if (front[39] < ulen && front[39] + 3 < vlen)
  {
    front[39] += suffix_or_prefix_match_len(
        useq, front[39], vseq, front[39] + 3, ulen, vlen, useqnum, vseqnum);
  }
  if (d_bound <= 3 && front[36 + d_bound] >= ulen)
  {
    return 3;
  }
  if (d_max == 3)
  {
    return 3 + 1;
  }
  /* d = 4: compute front[16...24] */
  /* h = -4 */
  front[0] = front[33] + 1;
  if (front[0] < ulen)
  {
    front[0] += suffix_or_prefix_match_len(useq, front[0], vseq, front[0] - 4,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -3 */
  front[1] = MAX2(front[33], front[34]) + 1;
  if (front[1] < ulen)
  {
    front[1] += suffix_or_prefix_match_len(useq, front[1], vseq, front[1] - 3,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -2 */
  front[2] = MAX2(front[33], MAX2(front[34], front[35]) + 1);
  if (front[2] < ulen)
  {
    front[2] += suffix_or_prefix_match_len(useq, front[2], vseq, front[2] - 2,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -1 */
  front[3] = MAX2(front[34], MAX2(front[35], front[36]) + 1);
  if (front[3] < ulen)
  {
    front[3] += suffix_or_prefix_match_len(useq, front[3], vseq, front[3] - 1,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 0 */
  front[4] = MAX2(front[35], MAX2(front[36], front[37]) + 1);
  if (front[4] < ulen)
  {
    front[4] += suffix_or_prefix_match_len(useq, front[4], vseq, front[4], ulen,
                                           vlen, useqnum, vseqnum);
  }
  /* h = 1 */
  front[5] = MAX2(front[36], MAX2(front[37], front[38]) + 1);
  if (front[5] < ulen && front[5] + 1 < vlen)
  {
    front[5] += suffix_or_prefix_match_len(useq, front[5], vseq, front[5] + 1,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 2 */
  front[6] = MAX2(front[37], MAX2(front[38], front[39]) + 1);
  if (front[6] < ulen && front[6] + 2 < vlen)
  {
    front[6] += suffix_or_prefix_match_len(useq, front[6], vseq, front[6] + 2,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 3 */
  front[7] = MAX2(front[38], front[39] + 1);
  if (front[7] < ulen && front[7] + 3 < vlen)
  {
    front[7] += suffix_or_prefix_match_len(useq, front[7], vseq, front[7] + 3,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 4 */
  front[8] = front[39];
  if (front[8] < ulen && front[8] + 4 < vlen)
  {
    front[8] += suffix_or_prefix_match_len(useq, front[8], vseq, front[8] + 4,
                                           ulen, vlen, useqnum, vseqnum);
  }
  if (d_bound <= 4 && front[4 + d_bound] >= ulen)
  {
    return 4;
  }
  if (d_max == 4)
  {
    return 4 + 1;
  }
  /* d = 5: compute front[25...35] */
  /* h = -5 */
  front[29] = front[0] + 1;
  if (front[29] < ulen)
  {
    front[29] += suffix_or_prefix_match_len(
        useq, front[29], vseq, front[29] - 5, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -4 */
  front[30] = MAX2(front[0], front[1]) + 1;
  if (front[30] < ulen)
  {
    front[30] += suffix_or_prefix_match_len(
        useq, front[30], vseq, front[30] - 4, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -3 */
  front[31] = MAX2(front[0], MAX2(front[1], front[2]) + 1);
  if (front[31] < ulen)
  {
    front[31] += suffix_or_prefix_match_len(
        useq, front[31], vseq, front[31] - 3, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -2 */
  front[32] = MAX2(front[1], MAX2(front[2], front[3]) + 1);
  if (front[32] < ulen)
  {
    front[32] += suffix_or_prefix_match_len(
        useq, front[32], vseq, front[32] - 2, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -1 */
  front[33] = MAX2(front[2], MAX2(front[3], front[4]) + 1);
  if (front[33] < ulen)
  {
    front[33] += suffix_or_prefix_match_len(
        useq, front[33], vseq, front[33] - 1, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 0 */
  front[34] = MAX2(front[3], MAX2(front[4], front[5]) + 1);
  if (front[34] < ulen)
  {
    front[34] += suffix_or_prefix_match_len(useq, front[34], vseq, front[34],
                                            ulen, vlen, useqnum, vseqnum);
  }
  /* h = 1 */
  front[35] = MAX2(front[4], MAX2(front[5], front[6]) + 1);
  if (front[35] < ulen && front[35] + 1 < vlen)
  {
    front[35] += suffix_or_prefix_match_len(
        useq, front[35], vseq, front[35] + 1, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 2 */
  front[36] = MAX2(front[5], MAX2(front[6], front[7]) + 1);
  if (front[36] < ulen && front[36] + 2 < vlen)
  {
    front[36] += suffix_or_prefix_match_len(
        useq, front[36], vseq, front[36] + 2, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 3 */
  front[37] = MAX2(front[6], MAX2(front[7], front[8]) + 1);
  if (front[37] < ulen && front[37] + 3 < vlen)
  {
    front[37] += suffix_or_prefix_match_len(
        useq, front[37], vseq, front[37] + 3, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 4 */
  front[38] = MAX2(front[7], front[8] + 1);
  if (front[38] < ulen && front[38] + 4 < vlen)
  {
    front[38] += suffix_or_prefix_match_len(
        useq, front[38], vseq, front[38] + 4, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 5 */
  front[39] = front[8];
  if (front[39] < ulen && front[39] + 5 < vlen)
  {
    front[39] += suffix_or_prefix_match_len(
        useq, front[39], vseq, front[39] + 5, ulen, vlen, useqnum, vseqnum);
  }
  if (d_bound <= 5 && front[34 + d_bound] >= ulen)
  {
    return 5;
  }
  if (d_max == 5)
  {
    return 5 + 1;
  }
  /* d = 6: compute front[36...48] */
  /* h = -6 */
  front[0] = front[29] + 1;
  if (front[0] < ulen)
  {
    front[0] += suffix_or_prefix_match_len(useq, front[0], vseq, front[0] - 6,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -5 */
  front[1] = MAX2(front[29], front[30]) + 1;
  if (front[1] < ulen)
  {
    front[1] += suffix_or_prefix_match_len(useq, front[1], vseq, front[1] - 5,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -4 */
  front[2] = MAX2(front[29], MAX2(front[30], front[31]) + 1);
  if (front[2] < ulen)
  {
    front[2] += suffix_or_prefix_match_len(useq, front[2], vseq, front[2] - 4,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -3 */
  front[3] = MAX2(front[30], MAX2(front[31], front[32]) + 1);
  if (front[3] < ulen)
  {
    front[3] += suffix_or_prefix_match_len(useq, front[3], vseq, front[3] - 3,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -2 */
  front[4] = MAX2(front[31], MAX2(front[32], front[33]) + 1);
  if (front[4] < ulen)
  {
    front[4] += suffix_or_prefix_match_len(useq, front[4], vseq, front[4] - 2,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -1 */
  front[5] = MAX2(front[32], MAX2(front[33], front[34]) + 1);
  if (front[5] < ulen)
  {
    front[5] += suffix_or_prefix_match_len(useq, front[5], vseq, front[5] - 1,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 0 */
  front[6] = MAX2(front[33], MAX2(front[34], front[35]) + 1);
  if (front[6] < ulen)
  {
    front[6] += suffix_or_prefix_match_len(useq, front[6], vseq, front[6], ulen,
                                           vlen, useqnum, vseqnum);
  }
  /* h = 1 */
  front[7] = MAX2(front[34], MAX2(front[35], front[36]) + 1);
  if (front[7] < ulen && front[7] + 1 < vlen)
  {
    front[7] += suffix_or_prefix_match_len(useq, front[7], vseq, front[7] + 1,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 2 */
  front[8] = MAX2(front[35], MAX2(front[36], front[37]) + 1);
  if (front[8] < ulen && front[8] + 2 < vlen)
  {
    front[8] += suffix_or_prefix_match_len(useq, front[8], vseq, front[8] + 2,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 3 */
  front[9] = MAX2(front[36], MAX2(front[37], front[38]) + 1);
  if (front[9] < ulen && front[9] + 3 < vlen)
  {
    front[9] += suffix_or_prefix_match_len(useq, front[9], vseq, front[9] + 3,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 4 */
  front[10] = MAX2(front[37], MAX2(front[38], front[39]) + 1);
  if (front[10] < ulen && front[10] + 4 < vlen)
  {
    front[10] += suffix_or_prefix_match_len(
        useq, front[10], vseq, front[10] + 4, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 5 */
  front[11] = MAX2(front[38], front[39] + 1);
  if (front[11] < ulen && front[11] + 5 < vlen)
  {
    front[11] += suffix_or_prefix_match_len(
        useq, front[11], vseq, front[11] + 5, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 6 */
  front[12] = front[39];
  if (front[12] < ulen && front[12] + 6 < vlen)
  {
    front[12] += suffix_or_prefix_match_len(
        useq, front[12], vseq, front[12] + 6, ulen, vlen, useqnum, vseqnum);
  }
  if (d_bound <= 6 && front[6 + d_bound] >= ulen)
  {
    return 6;
  }
  if (d_max == 6)
  {
    return 6 + 1;
  }
  /* d = 7: compute front[49...63] */
  /* h = -7 */
  front[25] = front[0] + 1;
  if (front[25] < ulen)
  {
    front[25] += suffix_or_prefix_match_len(
        useq, front[25], vseq, front[25] - 7, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -6 */
  front[26] = MAX2(front[0], front[1]) + 1;
  if (front[26] < ulen)
  {
    front[26] += suffix_or_prefix_match_len(
        useq, front[26], vseq, front[26] - 6, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -5 */
  front[27] = MAX2(front[0], MAX2(front[1], front[2]) + 1);
  if (front[27] < ulen)
  {
    front[27] += suffix_or_prefix_match_len(
        useq, front[27], vseq, front[27] - 5, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -4 */
  front[28] = MAX2(front[1], MAX2(front[2], front[3]) + 1);
  if (front[28] < ulen)
  {
    front[28] += suffix_or_prefix_match_len(
        useq, front[28], vseq, front[28] - 4, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -3 */
  front[29] = MAX2(front[2], MAX2(front[3], front[4]) + 1);
  if (front[29] < ulen)
  {
    front[29] += suffix_or_prefix_match_len(
        useq, front[29], vseq, front[29] - 3, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -2 */
  front[30] = MAX2(front[3], MAX2(front[4], front[5]) + 1);
  if (front[30] < ulen)
  {
    front[30] += suffix_or_prefix_match_len(
        useq, front[30], vseq, front[30] - 2, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -1 */
  front[31] = MAX2(front[4], MAX2(front[5], front[6]) + 1);
  if (front[31] < ulen)
  {
    front[31] += suffix_or_prefix_match_len(
        useq, front[31], vseq, front[31] - 1, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 0 */
  front[32] = MAX2(front[5], MAX2(front[6], front[7]) + 1);
  if (front[32] < ulen)
  {
    front[32] += suffix_or_prefix_match_len(useq, front[32], vseq, front[32],
                                            ulen, vlen, useqnum, vseqnum);
  }
  /* h = 1 */
  front[33] = MAX2(front[6], MAX2(front[7], front[8]) + 1);
  if (front[33] < ulen && front[33] + 1 < vlen)
  {
    front[33] += suffix_or_prefix_match_len(
        useq, front[33], vseq, front[33] + 1, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 2 */
  front[34] = MAX2(front[7], MAX2(front[8], front[9]) + 1);
  if (front[34] < ulen && front[34] + 2 < vlen)
  {
    front[34] += suffix_or_prefix_match_len(
        useq, front[34], vseq, front[34] + 2, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 3 */
  front[35] = MAX2(front[8], MAX2(front[9], front[10]) + 1);
  if (front[35] < ulen && front[35] + 3 < vlen)
  {
    front[35] += suffix_or_prefix_match_len(
        useq, front[35], vseq, front[35] + 3, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 4 */
  front[36] = MAX2(front[9], MAX2(front[10], front[11]) + 1);
  if (front[36] < ulen && front[36] + 4 < vlen)
  {
    front[36] += suffix_or_prefix_match_len(
        useq, front[36], vseq, front[36] + 4, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 5 */
  front[37] = MAX2(front[10], MAX2(front[11], front[12]) + 1);
  if (front[37] < ulen && front[37] + 5 < vlen)
  {
    front[37] += suffix_or_prefix_match_len(
        useq, front[37], vseq, front[37] + 5, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 6 */
  front[38] = MAX2(front[11], front[12] + 1);
  if (front[38] < ulen && front[38] + 6 < vlen)
  {
    front[38] += suffix_or_prefix_match_len(
        useq, front[38], vseq, front[38] + 6, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 7 */
  front[39] = front[12];
  if (front[39] < ulen && front[39] + 7 < vlen)
  {
    front[39] += suffix_or_prefix_match_len(
        useq, front[39], vseq, front[39] + 7, ulen, vlen, useqnum, vseqnum);
  }
  if (d_bound <= 7 && front[32 + d_bound] >= ulen)
  {
    return 7;
  }
  if (d_max == 7)
  {
    return 7 + 1;
  }
  /* d = 8: compute front[64...80] */
  /* h = -8 */
  front[0] = front[25] + 1;
  if (front[0] < ulen)
  {
    front[0] += suffix_or_prefix_match_len(useq, front[0], vseq, front[0] - 8,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -7 */
  front[1] = MAX2(front[25], front[26]) + 1;
  if (front[1] < ulen)
  {
    front[1] += suffix_or_prefix_match_len(useq, front[1], vseq, front[1] - 7,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -6 */
  front[2] = MAX2(front[25], MAX2(front[26], front[27]) + 1);
  if (front[2] < ulen)
  {
    front[2] += suffix_or_prefix_match_len(useq, front[2], vseq, front[2] - 6,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -5 */
  front[3] = MAX2(front[26], MAX2(front[27], front[28]) + 1);
  if (front[3] < ulen)
  {
    front[3] += suffix_or_prefix_match_len(useq, front[3], vseq, front[3] - 5,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -4 */
  front[4] = MAX2(front[27], MAX2(front[28], front[29]) + 1);
  if (front[4] < ulen)
  {
    front[4] += suffix_or_prefix_match_len(useq, front[4], vseq, front[4] - 4,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -3 */
  front[5] = MAX2(front[28], MAX2(front[29], front[30]) + 1);
  if (front[5] < ulen)
  {
    front[5] += suffix_or_prefix_match_len(useq, front[5], vseq, front[5] - 3,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -2 */
  front[6] = MAX2(front[29], MAX2(front[30], front[31]) + 1);
  if (front[6] < ulen)
  {
    front[6] += suffix_or_prefix_match_len(useq, front[6], vseq, front[6] - 2,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -1 */
  front[7] = MAX2(front[30], MAX2(front[31], front[32]) + 1);
  if (front[7] < ulen)
  {
    front[7] += suffix_or_prefix_match_len(useq, front[7], vseq, front[7] - 1,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 0 */
  front[8] = MAX2(front[31], MAX2(front[32], front[33]) + 1);
  if (front[8] < ulen)
  {
    front[8] += suffix_or_prefix_match_len(useq, front[8], vseq, front[8], ulen,
                                           vlen, useqnum, vseqnum);
  }
  /* h = 1 */
  front[9] = MAX2(front[32], MAX2(front[33], front[34]) + 1);
  if (front[9] < ulen && front[9] + 1 < vlen)
  {
    front[9] += suffix_or_prefix_match_len(useq, front[9], vseq, front[9] + 1,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 2 */
  front[10] = MAX2(front[33], MAX2(front[34], front[35]) + 1);
  if (front[10] < ulen && front[10] + 2 < vlen)
  {
    front[10] += suffix_or_prefix_match_len(
        useq, front[10], vseq, front[10] + 2, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 3 */
  front[11] = MAX2(front[34], MAX2(front[35], front[36]) + 1);
  if (front[11] < ulen && front[11] + 3 < vlen)
  {
    front[11] += suffix_or_prefix_match_len(
        useq, front[11], vseq, front[11] + 3, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 4 */
  front[12] = MAX2(front[35], MAX2(front[36], front[37]) + 1);
  if (front[12] < ulen && front[12] + 4 < vlen)
  {
    front[12] += suffix_or_prefix_match_len(
        useq, front[12], vseq, front[12] + 4, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 5 */
  front[13] = MAX2(front[36], MAX2(front[37], front[38]) + 1);
  if (front[13] < ulen && front[13] + 5 < vlen)
  {
    front[13] += suffix_or_prefix_match_len(
        useq, front[13], vseq, front[13] + 5, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 6 */
  front[14] = MAX2(front[37], MAX2(front[38], front[39]) + 1);
  if (front[14] < ulen && front[14] + 6 < vlen)
  {
    front[14] += suffix_or_prefix_match_len(
        useq, front[14], vseq, front[14] + 6, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 7 */
  front[15] = MAX2(front[38], front[39] + 1);
  if (front[15] < ulen && front[15] + 7 < vlen)
  {
    front[15] += suffix_or_prefix_match_len(
        useq, front[15], vseq, front[15] + 7, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 8 */
  front[16] = front[39];
  if (front[16] < ulen && front[16] + 8 < vlen)
  {
    front[16] += suffix_or_prefix_match_len(
        useq, front[16], vseq, front[16] + 8, ulen, vlen, useqnum, vseqnum);
  }
  if (d_bound <= 8 && front[8 + d_bound] >= ulen)
  {
    return 8;
  }
  if (d_max == 8)
  {
    return 8 + 1;
  }
  /* d = 9: compute front[81...99] */
  /* h = -9 */
  front[21] = front[0] + 1;
  if (front[21] < ulen)
  {
    front[21] += suffix_or_prefix_match_len(
        useq, front[21], vseq, front[21] - 9, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -8 */
  front[22] = MAX2(front[0], front[1]) + 1;
  if (front[22] < ulen)
  {
    front[22] += suffix_or_prefix_match_len(
        useq, front[22], vseq, front[22] - 8, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -7 */
  front[23] = MAX2(front[0], MAX2(front[1], front[2]) + 1);
  if (front[23] < ulen)
  {
    front[23] += suffix_or_prefix_match_len(
        useq, front[23], vseq, front[23] - 7, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -6 */
  front[24] = MAX2(front[1], MAX2(front[2], front[3]) + 1);
  if (front[24] < ulen)
  {
    front[24] += suffix_or_prefix_match_len(
        useq, front[24], vseq, front[24] - 6, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -5 */
  front[25] = MAX2(front[2], MAX2(front[3], front[4]) + 1);
  if (front[25] < ulen)
  {
    front[25] += suffix_or_prefix_match_len(
        useq, front[25], vseq, front[25] - 5, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -4 */
  front[26] = MAX2(front[3], MAX2(front[4], front[5]) + 1);
  if (front[26] < ulen)
  {
    front[26] += suffix_or_prefix_match_len(
        useq, front[26], vseq, front[26] - 4, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -3 */
  front[27] = MAX2(front[4], MAX2(front[5], front[6]) + 1);
  if (front[27] < ulen)
  {
    front[27] += suffix_or_prefix_match_len(
        useq, front[27], vseq, front[27] - 3, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -2 */
  front[28] = MAX2(front[5], MAX2(front[6], front[7]) + 1);
  if (front[28] < ulen)
  {
    front[28] += suffix_or_prefix_match_len(
        useq, front[28], vseq, front[28] - 2, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -1 */
  front[29] = MAX2(front[6], MAX2(front[7], front[8]) + 1);
  if (front[29] < ulen)
  {
    front[29] += suffix_or_prefix_match_len(
        useq, front[29], vseq, front[29] - 1, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 0 */
  front[30] = MAX2(front[7], MAX2(front[8], front[9]) + 1);
  if (front[30] < ulen)
  {
    front[30] += suffix_or_prefix_match_len(useq, front[30], vseq, front[30],
                                            ulen, vlen, useqnum, vseqnum);
  }
  /* h = 1 */
  front[31] = MAX2(front[8], MAX2(front[9], front[10]) + 1);
  if (front[31] < ulen && front[31] + 1 < vlen)
  {
    front[31] += suffix_or_prefix_match_len(
        useq, front[31], vseq, front[31] + 1, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 2 */
  front[32] = MAX2(front[9], MAX2(front[10], front[11]) + 1);
  if (front[32] < ulen && front[32] + 2 < vlen)
  {
    front[32] += suffix_or_prefix_match_len(
        useq, front[32], vseq, front[32] + 2, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 3 */
  front[33] = MAX2(front[10], MAX2(front[11], front[12]) + 1);
  if (front[33] < ulen && front[33] + 3 < vlen)
  {
    front[33] += suffix_or_prefix_match_len(
        useq, front[33], vseq, front[33] + 3, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 4 */
  front[34] = MAX2(front[11], MAX2(front[12], front[13]) + 1);
  if (front[34] < ulen && front[34] + 4 < vlen)
  {
    front[34] += suffix_or_prefix_match_len(
        useq, front[34], vseq, front[34] + 4, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 5 */
  front[35] = MAX2(front[12], MAX2(front[13], front[14]) + 1);
  if (front[35] < ulen && front[35] + 5 < vlen)
  {
    front[35] += suffix_or_prefix_match_len(
        useq, front[35], vseq, front[35] + 5, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 6 */
  front[36] = MAX2(front[13], MAX2(front[14], front[15]) + 1);
  if (front[36] < ulen && front[36] + 6 < vlen)
  {
    front[36] += suffix_or_prefix_match_len(
        useq, front[36], vseq, front[36] + 6, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 7 */
  front[37] = MAX2(front[14], MAX2(front[15], front[16]) + 1);
  if (front[37] < ulen && front[37] + 7 < vlen)
  {
    front[37] += suffix_or_prefix_match_len(
        useq, front[37], vseq, front[37] + 7, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 8 */
  front[38] = MAX2(front[15], front[16] + 1);
  if (front[38] < ulen && front[38] + 8 < vlen)
  {
    front[38] += suffix_or_prefix_match_len(
        useq, front[38], vseq, front[38] + 8, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 9 */
  front[39] = front[16];
  if (front[39] < ulen && front[39] + 9 < vlen)
  {
    front[39] += suffix_or_prefix_match_len(
        useq, front[39], vseq, front[39] + 9, ulen, vlen, useqnum, vseqnum);
  }
  if (d_bound <= 9 && front[30 + d_bound] >= ulen)
  {
    return 9;
  }
  if (d_max == 9)
  {
    return 9 + 1;
  }
  /* d = 10: compute front[100...120] */
  /* h = -10 */
  front[0] = front[21] + 1;
  if (front[0] < ulen)
  {
    front[0] += suffix_or_prefix_match_len(useq, front[0], vseq, front[0] - 10,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -9 */
  front[1] = MAX2(front[21], front[22]) + 1;
  if (front[1] < ulen)
  {
    front[1] += suffix_or_prefix_match_len(useq, front[1], vseq, front[1] - 9,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -8 */
  front[2] = MAX2(front[21], MAX2(front[22], front[23]) + 1);
  if (front[2] < ulen)
  {
    front[2] += suffix_or_prefix_match_len(useq, front[2], vseq, front[2] - 8,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -7 */
  front[3] = MAX2(front[22], MAX2(front[23], front[24]) + 1);
  if (front[3] < ulen)
  {
    front[3] += suffix_or_prefix_match_len(useq, front[3], vseq, front[3] - 7,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -6 */
  front[4] = MAX2(front[23], MAX2(front[24], front[25]) + 1);
  if (front[4] < ulen)
  {
    front[4] += suffix_or_prefix_match_len(useq, front[4], vseq, front[4] - 6,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -5 */
  front[5] = MAX2(front[24], MAX2(front[25], front[26]) + 1);
  if (front[5] < ulen)
  {
    front[5] += suffix_or_prefix_match_len(useq, front[5], vseq, front[5] - 5,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -4 */
  front[6] = MAX2(front[25], MAX2(front[26], front[27]) + 1);
  if (front[6] < ulen)
  {
    front[6] += suffix_or_prefix_match_len(useq, front[6], vseq, front[6] - 4,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -3 */
  front[7] = MAX2(front[26], MAX2(front[27], front[28]) + 1);
  if (front[7] < ulen)
  {
    front[7] += suffix_or_prefix_match_len(useq, front[7], vseq, front[7] - 3,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -2 */
  front[8] = MAX2(front[27], MAX2(front[28], front[29]) + 1);
  if (front[8] < ulen)
  {
    front[8] += suffix_or_prefix_match_len(useq, front[8], vseq, front[8] - 2,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -1 */
  front[9] = MAX2(front[28], MAX2(front[29], front[30]) + 1);
  if (front[9] < ulen)
  {
    front[9] += suffix_or_prefix_match_len(useq, front[9], vseq, front[9] - 1,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 0 */
  front[10] = MAX2(front[29], MAX2(front[30], front[31]) + 1);
  if (front[10] < ulen)
  {
    front[10] += suffix_or_prefix_match_len(useq, front[10], vseq, front[10],
                                            ulen, vlen, useqnum, vseqnum);
  }
  /* h = 1 */
  front[11] = MAX2(front[30], MAX2(front[31], front[32]) + 1);
  if (front[11] < ulen && front[11] + 1 < vlen)
  {
    front[11] += suffix_or_prefix_match_len(
        useq, front[11], vseq, front[11] + 1, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 2 */
  front[12] = MAX2(front[31], MAX2(front[32], front[33]) + 1);
  if (front[12] < ulen && front[12] + 2 < vlen)
  {
    front[12] += suffix_or_prefix_match_len(
        useq, front[12], vseq, front[12] + 2, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 3 */
  front[13] = MAX2(front[32], MAX2(front[33], front[34]) + 1);
  if (front[13] < ulen && front[13] + 3 < vlen)
  {
    front[13] += suffix_or_prefix_match_len(
        useq, front[13], vseq, front[13] + 3, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 4 */
  front[14] = MAX2(front[33], MAX2(front[34], front[35]) + 1);
  if (front[14] < ulen && front[14] + 4 < vlen)
  {
    front[14] += suffix_or_prefix_match_len(
        useq, front[14], vseq, front[14] + 4, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 5 */
  front[15] = MAX2(front[34], MAX2(front[35], front[36]) + 1);
  if (front[15] < ulen && front[15] + 5 < vlen)
  {
    front[15] += suffix_or_prefix_match_len(
        useq, front[15], vseq, front[15] + 5, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 6 */
  front[16] = MAX2(front[35], MAX2(front[36], front[37]) + 1);
  if (front[16] < ulen && front[16] + 6 < vlen)
  {
    front[16] += suffix_or_prefix_match_len(
        useq, front[16], vseq, front[16] + 6, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 7 */
  front[17] = MAX2(front[36], MAX2(front[37], front[38]) + 1);
  if (front[17] < ulen && front[17] + 7 < vlen)
  {
    front[17] += suffix_or_prefix_match_len(
        useq, front[17], vseq, front[17] + 7, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 8 */
  front[18] = MAX2(front[37], MAX2(front[38], front[39]) + 1);
  if (front[18] < ulen && front[18] + 8 < vlen)
  {
    front[18] += suffix_or_prefix_match_len(
        useq, front[18], vseq, front[18] + 8, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 9 */
  front[19] = MAX2(front[38], front[39] + 1);
  if (front[19] < ulen && front[19] + 9 < vlen)
  {
    front[19] += suffix_or_prefix_match_len(
        useq, front[19], vseq, front[19] + 9, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 10 */
  front[20] = front[39];
  if (front[20] < ulen && front[20] + 10 < vlen)
  {
    front[20] += suffix_or_prefix_match_len(
        useq, front[20], vseq, front[20] + 10, ulen, vlen, useqnum, vseqnum);
  }
  if (d_bound <= 10 && front[10 + d_bound] >= ulen)
  {
    return 10;
  }
  if (d_max == 10)
  {
    return 10 + 1;
  }
  return fastedist_inplace_continue<char_type, FrontValue, 10, true,
                                    suffix_or_prefix_match_len>(
      front + 0, d_max, useq, ulen, vseq, vlen, useqnum, vseqnum);
}
template <typename char_type, typename FrontValue,
          LcsLcpLenType suffix_or_prefix_match_len>
static size_t fastedist_unrolled_same_seq_length(
    size_t d_max, const char_type *useq, size_t ulen, const char_type *vseq,
    __attribute__((unused)) size_t vlen, size_t useqnum, size_t vseqnum)
{
  if (d_max == 0)
  {
    fprintf(stderr, "%s: illegal value d_max=0: can only handle d_max > 0\n",
            __func__);
    exit(EXIT_FAILURE);
  }
  assert(ulen == vlen);
  FrontValue front[40];
  front[0] = suffix_or_prefix_match_len(useq, 0, vseq, 0, ulen, vlen, useqnum,
                                        vseqnum);
  /* d = 1: compute front[1...3] */
  /* h = -1 */
  front[37] = front[0] + 1;
  if (front[37] < ulen)
  {
    front[37] += suffix_or_prefix_match_len(
        useq, front[37], vseq, front[37] - 1, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 0 */
  front[38] = front[0] + 1;
  if (front[38] < ulen)
  {
    front[38] += suffix_or_prefix_match_len(useq, front[38], vseq, front[38],
                                            ulen, vlen, useqnum, vseqnum);
  }
  /* h = 1 */
  front[39] = front[0];
  if (front[39] + 1 < ulen)
  {
    front[39] += suffix_or_prefix_match_len(
        useq, front[39], vseq, front[39] + 1, ulen, vlen, useqnum, vseqnum);
  }
  if (front[38] >= ulen)
  {
    return 1;
  }
  if (d_max == 1)
  {
    return 1 + 1;
  }
  /* d = 2: compute front[4...8] */
  /* h = -2 */
  front[0] = front[37] + 1;
  if (front[0] < ulen)
  {
    front[0] += suffix_or_prefix_match_len(useq, front[0], vseq, front[0] - 2,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -1 */
  front[1] = MAX2(front[37], front[38]) + 1;
  if (front[1] < ulen)
  {
    front[1] += suffix_or_prefix_match_len(useq, front[1], vseq, front[1] - 1,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 0 */
  front[2] = MAX2(front[37], MAX2(front[38], front[39]) + 1);
  if (front[2] < ulen)
  {
    front[2] += suffix_or_prefix_match_len(useq, front[2], vseq, front[2], ulen,
                                           vlen, useqnum, vseqnum);
  }
  /* h = 1 */
  front[3] = MAX2(front[38], front[39] + 1);
  if (front[3] + 1 < ulen)
  {
    front[3] += suffix_or_prefix_match_len(useq, front[3], vseq, front[3] + 1,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 2 */
  front[4] = front[39];
  if (front[4] + 2 < ulen)
  {
    front[4] += suffix_or_prefix_match_len(useq, front[4], vseq, front[4] + 2,
                                           ulen, vlen, useqnum, vseqnum);
  }
  if (front[2] >= ulen)
  {
    return 2;
  }
  if (d_max == 2)
  {
    return 2 + 1;
  }
  /* d = 3: compute front[9...15] */
  /* h = -3 */
  front[33] = front[0] + 1;
  if (front[33] < ulen)
  {
    front[33] += suffix_or_prefix_match_len(
        useq, front[33], vseq, front[33] - 3, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -2 */
  front[34] = MAX2(front[0], front[1]) + 1;
  if (front[34] < ulen)
  {
    front[34] += suffix_or_prefix_match_len(
        useq, front[34], vseq, front[34] - 2, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -1 */
  front[35] = MAX2(front[0], MAX2(front[1], front[2]) + 1);
  if (front[35] < ulen)
  {
    front[35] += suffix_or_prefix_match_len(
        useq, front[35], vseq, front[35] - 1, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 0 */
  front[36] = MAX2(front[1], MAX2(front[2], front[3]) + 1);
  if (front[36] < ulen)
  {
    front[36] += suffix_or_prefix_match_len(useq, front[36], vseq, front[36],
                                            ulen, vlen, useqnum, vseqnum);
  }
  /* h = 1 */
  front[37] = MAX2(front[2], MAX2(front[3], front[4]) + 1);
  if (front[37] + 1 < ulen)
  {
    front[37] += suffix_or_prefix_match_len(
        useq, front[37], vseq, front[37] + 1, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 2 */
  front[38] = MAX2(front[3], front[4] + 1);
  if (front[38] + 2 < ulen)
  {
    front[38] += suffix_or_prefix_match_len(
        useq, front[38], vseq, front[38] + 2, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 3 */
  front[39] = front[4];
  if (front[39] + 3 < ulen)
  {
    front[39] += suffix_or_prefix_match_len(
        useq, front[39], vseq, front[39] + 3, ulen, vlen, useqnum, vseqnum);
  }
  if (front[36] >= ulen)
  {
    return 3;
  }
  if (d_max == 3)
  {
    return 3 + 1;
  }
  /* d = 4: compute front[16...24] */
  /* h = -4 */
  front[0] = front[33] + 1;
  if (front[0] < ulen)
  {
    front[0] += suffix_or_prefix_match_len(useq, front[0], vseq, front[0] - 4,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -3 */
  front[1] = MAX2(front[33], front[34]) + 1;
  if (front[1] < ulen)
  {
    front[1] += suffix_or_prefix_match_len(useq, front[1], vseq, front[1] - 3,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -2 */
  front[2] = MAX2(front[33], MAX2(front[34], front[35]) + 1);
  if (front[2] < ulen)
  {
    front[2] += suffix_or_prefix_match_len(useq, front[2], vseq, front[2] - 2,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -1 */
  front[3] = MAX2(front[34], MAX2(front[35], front[36]) + 1);
  if (front[3] < ulen)
  {
    front[3] += suffix_or_prefix_match_len(useq, front[3], vseq, front[3] - 1,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 0 */
  front[4] = MAX2(front[35], MAX2(front[36], front[37]) + 1);
  if (front[4] < ulen)
  {
    front[4] += suffix_or_prefix_match_len(useq, front[4], vseq, front[4], ulen,
                                           vlen, useqnum, vseqnum);
  }
  /* h = 1 */
  front[5] = MAX2(front[36], MAX2(front[37], front[38]) + 1);
  if (front[5] + 1 < ulen)
  {
    front[5] += suffix_or_prefix_match_len(useq, front[5], vseq, front[5] + 1,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 2 */
  front[6] = MAX2(front[37], MAX2(front[38], front[39]) + 1);
  if (front[6] + 2 < ulen)
  {
    front[6] += suffix_or_prefix_match_len(useq, front[6], vseq, front[6] + 2,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 3 */
  front[7] = MAX2(front[38], front[39] + 1);
  if (front[7] + 3 < ulen)
  {
    front[7] += suffix_or_prefix_match_len(useq, front[7], vseq, front[7] + 3,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 4 */
  front[8] = front[39];
  if (front[8] + 4 < ulen)
  {
    front[8] += suffix_or_prefix_match_len(useq, front[8], vseq, front[8] + 4,
                                           ulen, vlen, useqnum, vseqnum);
  }
  if (front[4] >= ulen)
  {
    return 4;
  }
  if (d_max == 4)
  {
    return 4 + 1;
  }
  /* d = 5: compute front[25...35] */
  /* h = -5 */
  front[29] = front[0] + 1;
  if (front[29] < ulen)
  {
    front[29] += suffix_or_prefix_match_len(
        useq, front[29], vseq, front[29] - 5, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -4 */
  front[30] = MAX2(front[0], front[1]) + 1;
  if (front[30] < ulen)
  {
    front[30] += suffix_or_prefix_match_len(
        useq, front[30], vseq, front[30] - 4, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -3 */
  front[31] = MAX2(front[0], MAX2(front[1], front[2]) + 1);
  if (front[31] < ulen)
  {
    front[31] += suffix_or_prefix_match_len(
        useq, front[31], vseq, front[31] - 3, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -2 */
  front[32] = MAX2(front[1], MAX2(front[2], front[3]) + 1);
  if (front[32] < ulen)
  {
    front[32] += suffix_or_prefix_match_len(
        useq, front[32], vseq, front[32] - 2, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -1 */
  front[33] = MAX2(front[2], MAX2(front[3], front[4]) + 1);
  if (front[33] < ulen)
  {
    front[33] += suffix_or_prefix_match_len(
        useq, front[33], vseq, front[33] - 1, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 0 */
  front[34] = MAX2(front[3], MAX2(front[4], front[5]) + 1);
  if (front[34] < ulen)
  {
    front[34] += suffix_or_prefix_match_len(useq, front[34], vseq, front[34],
                                            ulen, vlen, useqnum, vseqnum);
  }
  /* h = 1 */
  front[35] = MAX2(front[4], MAX2(front[5], front[6]) + 1);
  if (front[35] + 1 < ulen)
  {
    front[35] += suffix_or_prefix_match_len(
        useq, front[35], vseq, front[35] + 1, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 2 */
  front[36] = MAX2(front[5], MAX2(front[6], front[7]) + 1);
  if (front[36] + 2 < ulen)
  {
    front[36] += suffix_or_prefix_match_len(
        useq, front[36], vseq, front[36] + 2, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 3 */
  front[37] = MAX2(front[6], MAX2(front[7], front[8]) + 1);
  if (front[37] + 3 < ulen)
  {
    front[37] += suffix_or_prefix_match_len(
        useq, front[37], vseq, front[37] + 3, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 4 */
  front[38] = MAX2(front[7], front[8] + 1);
  if (front[38] + 4 < ulen)
  {
    front[38] += suffix_or_prefix_match_len(
        useq, front[38], vseq, front[38] + 4, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 5 */
  front[39] = front[8];
  if (front[39] + 5 < ulen)
  {
    front[39] += suffix_or_prefix_match_len(
        useq, front[39], vseq, front[39] + 5, ulen, vlen, useqnum, vseqnum);
  }
  if (front[34] >= ulen)
  {
    return 5;
  }
  if (d_max == 5)
  {
    return 5 + 1;
  }
  /* d = 6: compute front[36...48] */
  /* h = -6 */
  front[0] = front[29] + 1;
  if (front[0] < ulen)
  {
    front[0] += suffix_or_prefix_match_len(useq, front[0], vseq, front[0] - 6,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -5 */
  front[1] = MAX2(front[29], front[30]) + 1;
  if (front[1] < ulen)
  {
    front[1] += suffix_or_prefix_match_len(useq, front[1], vseq, front[1] - 5,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -4 */
  front[2] = MAX2(front[29], MAX2(front[30], front[31]) + 1);
  if (front[2] < ulen)
  {
    front[2] += suffix_or_prefix_match_len(useq, front[2], vseq, front[2] - 4,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -3 */
  front[3] = MAX2(front[30], MAX2(front[31], front[32]) + 1);
  if (front[3] < ulen)
  {
    front[3] += suffix_or_prefix_match_len(useq, front[3], vseq, front[3] - 3,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -2 */
  front[4] = MAX2(front[31], MAX2(front[32], front[33]) + 1);
  if (front[4] < ulen)
  {
    front[4] += suffix_or_prefix_match_len(useq, front[4], vseq, front[4] - 2,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -1 */
  front[5] = MAX2(front[32], MAX2(front[33], front[34]) + 1);
  if (front[5] < ulen)
  {
    front[5] += suffix_or_prefix_match_len(useq, front[5], vseq, front[5] - 1,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 0 */
  front[6] = MAX2(front[33], MAX2(front[34], front[35]) + 1);
  if (front[6] < ulen)
  {
    front[6] += suffix_or_prefix_match_len(useq, front[6], vseq, front[6], ulen,
                                           vlen, useqnum, vseqnum);
  }
  /* h = 1 */
  front[7] = MAX2(front[34], MAX2(front[35], front[36]) + 1);
  if (front[7] + 1 < ulen)
  {
    front[7] += suffix_or_prefix_match_len(useq, front[7], vseq, front[7] + 1,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 2 */
  front[8] = MAX2(front[35], MAX2(front[36], front[37]) + 1);
  if (front[8] + 2 < ulen)
  {
    front[8] += suffix_or_prefix_match_len(useq, front[8], vseq, front[8] + 2,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 3 */
  front[9] = MAX2(front[36], MAX2(front[37], front[38]) + 1);
  if (front[9] + 3 < ulen)
  {
    front[9] += suffix_or_prefix_match_len(useq, front[9], vseq, front[9] + 3,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 4 */
  front[10] = MAX2(front[37], MAX2(front[38], front[39]) + 1);
  if (front[10] + 4 < ulen)
  {
    front[10] += suffix_or_prefix_match_len(
        useq, front[10], vseq, front[10] + 4, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 5 */
  front[11] = MAX2(front[38], front[39] + 1);
  if (front[11] + 5 < ulen)
  {
    front[11] += suffix_or_prefix_match_len(
        useq, front[11], vseq, front[11] + 5, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 6 */
  front[12] = front[39];
  if (front[12] + 6 < ulen)
  {
    front[12] += suffix_or_prefix_match_len(
        useq, front[12], vseq, front[12] + 6, ulen, vlen, useqnum, vseqnum);
  }
  if (front[6] >= ulen)
  {
    return 6;
  }
  if (d_max == 6)
  {
    return 6 + 1;
  }
  /* d = 7: compute front[49...63] */
  /* h = -7 */
  front[25] = front[0] + 1;
  if (front[25] < ulen)
  {
    front[25] += suffix_or_prefix_match_len(
        useq, front[25], vseq, front[25] - 7, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -6 */
  front[26] = MAX2(front[0], front[1]) + 1;
  if (front[26] < ulen)
  {
    front[26] += suffix_or_prefix_match_len(
        useq, front[26], vseq, front[26] - 6, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -5 */
  front[27] = MAX2(front[0], MAX2(front[1], front[2]) + 1);
  if (front[27] < ulen)
  {
    front[27] += suffix_or_prefix_match_len(
        useq, front[27], vseq, front[27] - 5, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -4 */
  front[28] = MAX2(front[1], MAX2(front[2], front[3]) + 1);
  if (front[28] < ulen)
  {
    front[28] += suffix_or_prefix_match_len(
        useq, front[28], vseq, front[28] - 4, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -3 */
  front[29] = MAX2(front[2], MAX2(front[3], front[4]) + 1);
  if (front[29] < ulen)
  {
    front[29] += suffix_or_prefix_match_len(
        useq, front[29], vseq, front[29] - 3, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -2 */
  front[30] = MAX2(front[3], MAX2(front[4], front[5]) + 1);
  if (front[30] < ulen)
  {
    front[30] += suffix_or_prefix_match_len(
        useq, front[30], vseq, front[30] - 2, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -1 */
  front[31] = MAX2(front[4], MAX2(front[5], front[6]) + 1);
  if (front[31] < ulen)
  {
    front[31] += suffix_or_prefix_match_len(
        useq, front[31], vseq, front[31] - 1, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 0 */
  front[32] = MAX2(front[5], MAX2(front[6], front[7]) + 1);
  if (front[32] < ulen)
  {
    front[32] += suffix_or_prefix_match_len(useq, front[32], vseq, front[32],
                                            ulen, vlen, useqnum, vseqnum);
  }
  /* h = 1 */
  front[33] = MAX2(front[6], MAX2(front[7], front[8]) + 1);
  if (front[33] + 1 < ulen)
  {
    front[33] += suffix_or_prefix_match_len(
        useq, front[33], vseq, front[33] + 1, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 2 */
  front[34] = MAX2(front[7], MAX2(front[8], front[9]) + 1);
  if (front[34] + 2 < ulen)
  {
    front[34] += suffix_or_prefix_match_len(
        useq, front[34], vseq, front[34] + 2, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 3 */
  front[35] = MAX2(front[8], MAX2(front[9], front[10]) + 1);
  if (front[35] + 3 < ulen)
  {
    front[35] += suffix_or_prefix_match_len(
        useq, front[35], vseq, front[35] + 3, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 4 */
  front[36] = MAX2(front[9], MAX2(front[10], front[11]) + 1);
  if (front[36] + 4 < ulen)
  {
    front[36] += suffix_or_prefix_match_len(
        useq, front[36], vseq, front[36] + 4, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 5 */
  front[37] = MAX2(front[10], MAX2(front[11], front[12]) + 1);
  if (front[37] + 5 < ulen)
  {
    front[37] += suffix_or_prefix_match_len(
        useq, front[37], vseq, front[37] + 5, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 6 */
  front[38] = MAX2(front[11], front[12] + 1);
  if (front[38] + 6 < ulen)
  {
    front[38] += suffix_or_prefix_match_len(
        useq, front[38], vseq, front[38] + 6, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 7 */
  front[39] = front[12];
  if (front[39] + 7 < ulen)
  {
    front[39] += suffix_or_prefix_match_len(
        useq, front[39], vseq, front[39] + 7, ulen, vlen, useqnum, vseqnum);
  }
  if (front[32] >= ulen)
  {
    return 7;
  }
  if (d_max == 7)
  {
    return 7 + 1;
  }
  /* d = 8: compute front[64...80] */
  /* h = -8 */
  front[0] = front[25] + 1;
  if (front[0] < ulen)
  {
    front[0] += suffix_or_prefix_match_len(useq, front[0], vseq, front[0] - 8,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -7 */
  front[1] = MAX2(front[25], front[26]) + 1;
  if (front[1] < ulen)
  {
    front[1] += suffix_or_prefix_match_len(useq, front[1], vseq, front[1] - 7,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -6 */
  front[2] = MAX2(front[25], MAX2(front[26], front[27]) + 1);
  if (front[2] < ulen)
  {
    front[2] += suffix_or_prefix_match_len(useq, front[2], vseq, front[2] - 6,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -5 */
  front[3] = MAX2(front[26], MAX2(front[27], front[28]) + 1);
  if (front[3] < ulen)
  {
    front[3] += suffix_or_prefix_match_len(useq, front[3], vseq, front[3] - 5,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -4 */
  front[4] = MAX2(front[27], MAX2(front[28], front[29]) + 1);
  if (front[4] < ulen)
  {
    front[4] += suffix_or_prefix_match_len(useq, front[4], vseq, front[4] - 4,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -3 */
  front[5] = MAX2(front[28], MAX2(front[29], front[30]) + 1);
  if (front[5] < ulen)
  {
    front[5] += suffix_or_prefix_match_len(useq, front[5], vseq, front[5] - 3,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -2 */
  front[6] = MAX2(front[29], MAX2(front[30], front[31]) + 1);
  if (front[6] < ulen)
  {
    front[6] += suffix_or_prefix_match_len(useq, front[6], vseq, front[6] - 2,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -1 */
  front[7] = MAX2(front[30], MAX2(front[31], front[32]) + 1);
  if (front[7] < ulen)
  {
    front[7] += suffix_or_prefix_match_len(useq, front[7], vseq, front[7] - 1,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 0 */
  front[8] = MAX2(front[31], MAX2(front[32], front[33]) + 1);
  if (front[8] < ulen)
  {
    front[8] += suffix_or_prefix_match_len(useq, front[8], vseq, front[8], ulen,
                                           vlen, useqnum, vseqnum);
  }
  /* h = 1 */
  front[9] = MAX2(front[32], MAX2(front[33], front[34]) + 1);
  if (front[9] + 1 < ulen)
  {
    front[9] += suffix_or_prefix_match_len(useq, front[9], vseq, front[9] + 1,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 2 */
  front[10] = MAX2(front[33], MAX2(front[34], front[35]) + 1);
  if (front[10] + 2 < ulen)
  {
    front[10] += suffix_or_prefix_match_len(
        useq, front[10], vseq, front[10] + 2, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 3 */
  front[11] = MAX2(front[34], MAX2(front[35], front[36]) + 1);
  if (front[11] + 3 < ulen)
  {
    front[11] += suffix_or_prefix_match_len(
        useq, front[11], vseq, front[11] + 3, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 4 */
  front[12] = MAX2(front[35], MAX2(front[36], front[37]) + 1);
  if (front[12] + 4 < ulen)
  {
    front[12] += suffix_or_prefix_match_len(
        useq, front[12], vseq, front[12] + 4, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 5 */
  front[13] = MAX2(front[36], MAX2(front[37], front[38]) + 1);
  if (front[13] + 5 < ulen)
  {
    front[13] += suffix_or_prefix_match_len(
        useq, front[13], vseq, front[13] + 5, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 6 */
  front[14] = MAX2(front[37], MAX2(front[38], front[39]) + 1);
  if (front[14] + 6 < ulen)
  {
    front[14] += suffix_or_prefix_match_len(
        useq, front[14], vseq, front[14] + 6, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 7 */
  front[15] = MAX2(front[38], front[39] + 1);
  if (front[15] + 7 < ulen)
  {
    front[15] += suffix_or_prefix_match_len(
        useq, front[15], vseq, front[15] + 7, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 8 */
  front[16] = front[39];
  if (front[16] + 8 < ulen)
  {
    front[16] += suffix_or_prefix_match_len(
        useq, front[16], vseq, front[16] + 8, ulen, vlen, useqnum, vseqnum);
  }
  if (front[8] >= ulen)
  {
    return 8;
  }
  if (d_max == 8)
  {
    return 8 + 1;
  }
  /* d = 9: compute front[81...99] */
  /* h = -9 */
  front[21] = front[0] + 1;
  if (front[21] < ulen)
  {
    front[21] += suffix_or_prefix_match_len(
        useq, front[21], vseq, front[21] - 9, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -8 */
  front[22] = MAX2(front[0], front[1]) + 1;
  if (front[22] < ulen)
  {
    front[22] += suffix_or_prefix_match_len(
        useq, front[22], vseq, front[22] - 8, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -7 */
  front[23] = MAX2(front[0], MAX2(front[1], front[2]) + 1);
  if (front[23] < ulen)
  {
    front[23] += suffix_or_prefix_match_len(
        useq, front[23], vseq, front[23] - 7, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -6 */
  front[24] = MAX2(front[1], MAX2(front[2], front[3]) + 1);
  if (front[24] < ulen)
  {
    front[24] += suffix_or_prefix_match_len(
        useq, front[24], vseq, front[24] - 6, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -5 */
  front[25] = MAX2(front[2], MAX2(front[3], front[4]) + 1);
  if (front[25] < ulen)
  {
    front[25] += suffix_or_prefix_match_len(
        useq, front[25], vseq, front[25] - 5, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -4 */
  front[26] = MAX2(front[3], MAX2(front[4], front[5]) + 1);
  if (front[26] < ulen)
  {
    front[26] += suffix_or_prefix_match_len(
        useq, front[26], vseq, front[26] - 4, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -3 */
  front[27] = MAX2(front[4], MAX2(front[5], front[6]) + 1);
  if (front[27] < ulen)
  {
    front[27] += suffix_or_prefix_match_len(
        useq, front[27], vseq, front[27] - 3, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -2 */
  front[28] = MAX2(front[5], MAX2(front[6], front[7]) + 1);
  if (front[28] < ulen)
  {
    front[28] += suffix_or_prefix_match_len(
        useq, front[28], vseq, front[28] - 2, ulen, vlen, useqnum, vseqnum);
  }
  /* h = -1 */
  front[29] = MAX2(front[6], MAX2(front[7], front[8]) + 1);
  if (front[29] < ulen)
  {
    front[29] += suffix_or_prefix_match_len(
        useq, front[29], vseq, front[29] - 1, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 0 */
  front[30] = MAX2(front[7], MAX2(front[8], front[9]) + 1);
  if (front[30] < ulen)
  {
    front[30] += suffix_or_prefix_match_len(useq, front[30], vseq, front[30],
                                            ulen, vlen, useqnum, vseqnum);
  }
  /* h = 1 */
  front[31] = MAX2(front[8], MAX2(front[9], front[10]) + 1);
  if (front[31] + 1 < ulen)
  {
    front[31] += suffix_or_prefix_match_len(
        useq, front[31], vseq, front[31] + 1, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 2 */
  front[32] = MAX2(front[9], MAX2(front[10], front[11]) + 1);
  if (front[32] + 2 < ulen)
  {
    front[32] += suffix_or_prefix_match_len(
        useq, front[32], vseq, front[32] + 2, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 3 */
  front[33] = MAX2(front[10], MAX2(front[11], front[12]) + 1);
  if (front[33] + 3 < ulen)
  {
    front[33] += suffix_or_prefix_match_len(
        useq, front[33], vseq, front[33] + 3, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 4 */
  front[34] = MAX2(front[11], MAX2(front[12], front[13]) + 1);
  if (front[34] + 4 < ulen)
  {
    front[34] += suffix_or_prefix_match_len(
        useq, front[34], vseq, front[34] + 4, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 5 */
  front[35] = MAX2(front[12], MAX2(front[13], front[14]) + 1);
  if (front[35] + 5 < ulen)
  {
    front[35] += suffix_or_prefix_match_len(
        useq, front[35], vseq, front[35] + 5, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 6 */
  front[36] = MAX2(front[13], MAX2(front[14], front[15]) + 1);
  if (front[36] + 6 < ulen)
  {
    front[36] += suffix_or_prefix_match_len(
        useq, front[36], vseq, front[36] + 6, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 7 */
  front[37] = MAX2(front[14], MAX2(front[15], front[16]) + 1);
  if (front[37] + 7 < ulen)
  {
    front[37] += suffix_or_prefix_match_len(
        useq, front[37], vseq, front[37] + 7, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 8 */
  front[38] = MAX2(front[15], front[16] + 1);
  if (front[38] + 8 < ulen)
  {
    front[38] += suffix_or_prefix_match_len(
        useq, front[38], vseq, front[38] + 8, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 9 */
  front[39] = front[16];
  if (front[39] + 9 < ulen)
  {
    front[39] += suffix_or_prefix_match_len(
        useq, front[39], vseq, front[39] + 9, ulen, vlen, useqnum, vseqnum);
  }
  if (front[30] >= ulen)
  {
    return 9;
  }
  if (d_max == 9)
  {
    return 9 + 1;
  }
  /* d = 10: compute front[100...120] */
  /* h = -10 */
  front[0] = front[21] + 1;
  if (front[0] < ulen)
  {
    front[0] += suffix_or_prefix_match_len(useq, front[0], vseq, front[0] - 10,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -9 */
  front[1] = MAX2(front[21], front[22]) + 1;
  if (front[1] < ulen)
  {
    front[1] += suffix_or_prefix_match_len(useq, front[1], vseq, front[1] - 9,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -8 */
  front[2] = MAX2(front[21], MAX2(front[22], front[23]) + 1);
  if (front[2] < ulen)
  {
    front[2] += suffix_or_prefix_match_len(useq, front[2], vseq, front[2] - 8,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -7 */
  front[3] = MAX2(front[22], MAX2(front[23], front[24]) + 1);
  if (front[3] < ulen)
  {
    front[3] += suffix_or_prefix_match_len(useq, front[3], vseq, front[3] - 7,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -6 */
  front[4] = MAX2(front[23], MAX2(front[24], front[25]) + 1);
  if (front[4] < ulen)
  {
    front[4] += suffix_or_prefix_match_len(useq, front[4], vseq, front[4] - 6,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -5 */
  front[5] = MAX2(front[24], MAX2(front[25], front[26]) + 1);
  if (front[5] < ulen)
  {
    front[5] += suffix_or_prefix_match_len(useq, front[5], vseq, front[5] - 5,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -4 */
  front[6] = MAX2(front[25], MAX2(front[26], front[27]) + 1);
  if (front[6] < ulen)
  {
    front[6] += suffix_or_prefix_match_len(useq, front[6], vseq, front[6] - 4,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -3 */
  front[7] = MAX2(front[26], MAX2(front[27], front[28]) + 1);
  if (front[7] < ulen)
  {
    front[7] += suffix_or_prefix_match_len(useq, front[7], vseq, front[7] - 3,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -2 */
  front[8] = MAX2(front[27], MAX2(front[28], front[29]) + 1);
  if (front[8] < ulen)
  {
    front[8] += suffix_or_prefix_match_len(useq, front[8], vseq, front[8] - 2,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = -1 */
  front[9] = MAX2(front[28], MAX2(front[29], front[30]) + 1);
  if (front[9] < ulen)
  {
    front[9] += suffix_or_prefix_match_len(useq, front[9], vseq, front[9] - 1,
                                           ulen, vlen, useqnum, vseqnum);
  }
  /* h = 0 */
  front[10] = MAX2(front[29], MAX2(front[30], front[31]) + 1);
  if (front[10] < ulen)
  {
    front[10] += suffix_or_prefix_match_len(useq, front[10], vseq, front[10],
                                            ulen, vlen, useqnum, vseqnum);
  }
  /* h = 1 */
  front[11] = MAX2(front[30], MAX2(front[31], front[32]) + 1);
  if (front[11] + 1 < ulen)
  {
    front[11] += suffix_or_prefix_match_len(
        useq, front[11], vseq, front[11] + 1, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 2 */
  front[12] = MAX2(front[31], MAX2(front[32], front[33]) + 1);
  if (front[12] + 2 < ulen)
  {
    front[12] += suffix_or_prefix_match_len(
        useq, front[12], vseq, front[12] + 2, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 3 */
  front[13] = MAX2(front[32], MAX2(front[33], front[34]) + 1);
  if (front[13] + 3 < ulen)
  {
    front[13] += suffix_or_prefix_match_len(
        useq, front[13], vseq, front[13] + 3, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 4 */
  front[14] = MAX2(front[33], MAX2(front[34], front[35]) + 1);
  if (front[14] + 4 < ulen)
  {
    front[14] += suffix_or_prefix_match_len(
        useq, front[14], vseq, front[14] + 4, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 5 */
  front[15] = MAX2(front[34], MAX2(front[35], front[36]) + 1);
  if (front[15] + 5 < ulen)
  {
    front[15] += suffix_or_prefix_match_len(
        useq, front[15], vseq, front[15] + 5, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 6 */
  front[16] = MAX2(front[35], MAX2(front[36], front[37]) + 1);
  if (front[16] + 6 < ulen)
  {
    front[16] += suffix_or_prefix_match_len(
        useq, front[16], vseq, front[16] + 6, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 7 */
  front[17] = MAX2(front[36], MAX2(front[37], front[38]) + 1);
  if (front[17] + 7 < ulen)
  {
    front[17] += suffix_or_prefix_match_len(
        useq, front[17], vseq, front[17] + 7, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 8 */
  front[18] = MAX2(front[37], MAX2(front[38], front[39]) + 1);
  if (front[18] + 8 < ulen)
  {
    front[18] += suffix_or_prefix_match_len(
        useq, front[18], vseq, front[18] + 8, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 9 */
  front[19] = MAX2(front[38], front[39] + 1);
  if (front[19] + 9 < ulen)
  {
    front[19] += suffix_or_prefix_match_len(
        useq, front[19], vseq, front[19] + 9, ulen, vlen, useqnum, vseqnum);
  }
  /* h = 10 */
  front[20] = front[39];
  if (front[20] + 10 < ulen)
  {
    front[20] += suffix_or_prefix_match_len(
        useq, front[20], vseq, front[20] + 10, ulen, vlen, useqnum, vseqnum);
  }
  if (front[10] >= ulen)
  {
    return 10;
  }
  if (d_max == 10)
  {
    return 10 + 1;
  }
  return fastedist_inplace_continue<char_type, FrontValue, 10, true,
                                    suffix_or_prefix_match_len>(
      front + 0, d_max, useq, ulen, vseq, vlen, useqnum, vseqnum);
}
#endif
