/*
  Copyright (c) 2021 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2021 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/
/* This file is generated by ./split_num.py. Do not edit. */
#ifndef BITS2SPLIT_HPP
#define BITS2SPLIT_HPP
#include <cassert>
#include <cstddef>
#define MIN_NUM_SPLIT 1
#define MAX_NUM_SPLIT 64
static const int bits2split[] = {
  /* 1 */ 1,
  /* 2 */ 2,
  /* 3 */ 3,
  /* 4 */ 4,
  /* 5 */ 5,
  /* 6 */ 6,
  /* 7 */ 7,
  /* 8 */ 8,
  /* 9 */ 9,
  /* 10 */ 5, 5,
  /* 11 */ 5, 6,
  /* 12 */ 6, 6,
  /* 13 */ 6, 7,
  /* 14 */ 7, 7,
  /* 15 */ 7, 8,
  /* 16 */ 8, 8,
  /* 17 */ 8, 9,
  /* 18 */ 9, 9,
  /* 19 */ 6, 6, 7,
  /* 20 */ 6, 7, 7,
  /* 21 */ 7, 7, 7,
  /* 22 */ 7, 7, 8,
  /* 23 */ 7, 8, 8,
  /* 24 */ 8, 8, 8,
  /* 25 */ 8, 8, 9,
  /* 26 */ 8, 9, 9,
  /* 27 */ 9, 9, 9,
  /* 28 */ 7, 7, 7, 7,
  /* 29 */ 7, 7, 7, 8,
  /* 30 */ 7, 7, 8, 8,
  /* 31 */ 7, 8, 8, 8,
  /* 32 */ 8, 8, 8, 8,
  /* 33 */ 8, 8, 8, 9,
  /* 34 */ 8, 8, 9, 9,
  /* 35 */ 8, 9, 9, 9,
  /* 36 */ 9, 9, 9, 9,
  /* 37 */ 7, 7, 7, 8, 8,
  /* 38 */ 7, 7, 8, 8, 8,
  /* 39 */ 7, 8, 8, 8, 8,
  /* 40 */ 8, 8, 8, 8, 8,
  /* 41 */ 8, 8, 8, 8, 9,
  /* 42 */ 8, 8, 8, 9, 9,
  /* 43 */ 8, 8, 9, 9, 9,
  /* 44 */ 8, 9, 9, 9, 9,
  /* 45 */ 9, 9, 9, 9, 9,
  /* 46 */ 7, 7, 8, 8, 8, 8,
  /* 47 */ 7, 8, 8, 8, 8, 8,
  /* 48 */ 8, 8, 8, 8, 8, 8,
  /* 49 */ 8, 8, 8, 8, 8, 9,
  /* 50 */ 8, 8, 8, 8, 9, 9,
  /* 51 */ 8, 8, 8, 9, 9, 9,
  /* 52 */ 8, 8, 9, 9, 9, 9,
  /* 53 */ 8, 9, 9, 9, 9, 9,
  /* 54 */ 9, 9, 9, 9, 9, 9,
  /* 55 */ 7, 8, 8, 8, 8, 8, 8,
  /* 56 */ 8, 8, 8, 8, 8, 8, 8,
  /* 57 */ 8, 8, 8, 8, 8, 8, 9,
  /* 58 */ 8, 8, 8, 8, 8, 9, 9,
  /* 59 */ 8, 8, 8, 8, 9, 9, 9,
  /* 60 */ 8, 8, 8, 9, 9, 9, 9,
  /* 61 */ 8, 8, 9, 9, 9, 9, 9,
  /* 62 */ 8, 9, 9, 9, 9, 9, 9,
  /* 63 */ 9, 9, 9, 9, 9, 9, 9,
  /* 64 */ 8, 8, 8, 8, 8, 8, 8, 8
};
const int split_group_start[] = {
0, /* 1 */
1, /* 2 */
2, /* 3 */
3, /* 4 */
4, /* 5 */
5, /* 6 */
6, /* 7 */
7, /* 8 */
8, /* 9 */
9, /* 10 */
11, /* 11 */
13, /* 12 */
15, /* 13 */
17, /* 14 */
19, /* 15 */
21, /* 16 */
23, /* 17 */
25, /* 18 */
27, /* 19 */
30, /* 20 */
33, /* 21 */
36, /* 22 */
39, /* 23 */
42, /* 24 */
45, /* 25 */
48, /* 26 */
51, /* 27 */
54, /* 28 */
58, /* 29 */
62, /* 30 */
66, /* 31 */
70, /* 32 */
74, /* 33 */
78, /* 34 */
82, /* 35 */
86, /* 36 */
90, /* 37 */
95, /* 38 */
100, /* 39 */
105, /* 40 */
110, /* 41 */
115, /* 42 */
120, /* 43 */
125, /* 44 */
130, /* 45 */
135, /* 46 */
141, /* 47 */
147, /* 48 */
153, /* 49 */
159, /* 50 */
165, /* 51 */
171, /* 52 */
177, /* 53 */
183, /* 54 */
189, /* 55 */
196, /* 56 */
203, /* 57 */
210, /* 58 */
217, /* 59 */
224, /* 60 */
231, /* 61 */
238, /* 62 */
245, /* 63 */
252, /* 64 */
260
};
#define BIT_COUNTS_LAST_NUM_RANGE 8
static inline const int *bit_counts_get(size_t *num_ranges, int numbits)
{
#ifndef NDEBUG
  const size_t num_groups = sizeof split_group_start/
                            sizeof split_group_start[0],
               num_split_elems = sizeof bits2split/
                                 sizeof bits2split[0];
#endif
  assert(numbits >= MIN_NUM_SPLIT);
  assert(numbits <= MAX_NUM_SPLIT);
  assert(split_group_start[num_groups-2] + BIT_COUNTS_LAST_NUM_RANGE
           == num_split_elems);
  const int start = split_group_start[numbits - MIN_NUM_SPLIT];
  assert(static_cast<size_t>(numbits + 1 - MIN_NUM_SPLIT) < num_groups);
  const int end = split_group_start[numbits + 1 - MIN_NUM_SPLIT];
  *num_ranges = static_cast<size_t>(end - start);
  return bits2split + start;
}
#endif
