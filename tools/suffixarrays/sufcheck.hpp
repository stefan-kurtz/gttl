/*
  Copyright (c) 2012-2025 Stefan Kurtz <stefan.kurtz@uni-hamburg.de>
  Copyright (c) 2012-2025 Center for Bioinformatics, University of Hamburg

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

#ifndef SUFCHECK_HPP
#define SUFCHECK_HPP

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "utilities/multibitvector.hpp"

class Rangewithchar
{
  const size_t start;
#ifndef NDEBUG
  const size_t end;
#endif
  const uint8_t firstchar;
  public:
  Rangewithchar(size_t _start, [[maybe_unused]] size_t _end, uint8_t _firstchar)
    : start(_start)
#ifndef NDEBUG
    , end(_end)
#endif
    , firstchar(_firstchar)
  { }
  [[nodiscard]] size_t start_get(void) const { return start; }
#ifndef NDEBUG
  size_t width_get(void) const
  {
    return end - start + 1;
  }
#endif
  [[nodiscard]] uint8_t firstchar_get(void) const { return firstchar; }
};

template <bool is_multiseq>
static inline bool isspecial(uint8_t padding_char, uint8_t cc, uint8_t wildcard)
{
  static_assert(is_multiseq);
  return cc == padding_char or cc == wildcard;
}

/* The following funktion implements the linear time algorithm of
   @INPROCEEDINGS{BUR:KAER:2003,
   author = {Burkhardt, S. and K{\"a}rkk{\"a}inen, J.},
   title = {{Fast Lightweight Suffix Array Construction and Checking}},
   booktitle = {{Proceedings of the 14th Annual Symposium on Combinatorial
                 Pattern Matching (CPM)}},
   year = {2003},
   editor = {{Baeza-Yates, R. and Ch{\'a}vez, E. and Crochemore, M.}},
   volume = {2676},
   series = {LNCS},
   pages = {200-210},
   publisher = {Springer-Verlag}
   }
   to check the following suffix-order condition of the sorted suffix array:
   For all characters c, if SA[i,j] contains the suffixes starting
   with charcter c, then SA[i]+1, SA[i+1]+1, \ldots, SA[j]+1 occur
   in SA in this order (but not consecutively in general).
   The running time of the algorithm is independent of the alphabet size.
   The main problem is that it requires random access to the sequence
   which slows it down.
*/
template <typename SuftabBaseType, typename SequenceType, bool is_multiseq>
static void gttl_suftab_bk_suffixorder(const SequenceType *sequence,
                                       size_t totallength,
                                       [[maybe_unused]] uint8_t wildcard,
                                       [[maybe_unused]] uint8_t padding_char,
                                       size_t numofchars,
                                       const SuftabBaseType *suftab,
                                       const std::vector<Rangewithchar>
                                         &rangestore)
{
  size_t *const nexttab = static_cast<size_t *>(
                               calloc(numofchars, sizeof *nexttab));

  for (auto && rng : rangestore)
  {
    nexttab[rng.firstchar_get()] = rng.start_get();
  }
  for (size_t idx = 0; idx < totallength; idx++)
  {
    const SuftabBaseType position = suftab[idx];

    if (position > 0)
    {
      assert(position - 1 < totallength);

      if constexpr (is_multiseq)
      {
        const uint8_t cc = static_cast<uint8_t>(sequence[position - 1]);
        if (!isspecial<is_multiseq>(padding_char, cc, wildcard))
        {
          const SuftabBaseType checkpos
            = suftab[nexttab[static_cast<int>(cc)]] + 1;
          if (checkpos != position)
          {
            fprintf(stderr, "idx=%zu,checkpos=%zu,position=%zu\n", idx,
                    static_cast<size_t>(checkpos),
                    static_cast<size_t>(position));
            exit(EXIT_FAILURE);
          }
          nexttab[static_cast<int>(cc)]++;
        }
      } else
      {
        const uint8_t cc = static_cast<uint8_t>(sequence[position - 1]);
        const SuftabBaseType checkpos
          = suftab[nexttab[static_cast<int>(cc)]] + 1;
        if (checkpos != position)
        {
          fprintf(stderr, "idx=%zu,checkpos=%zu,position=%zu\n", idx,
                  static_cast<size_t>(checkpos),
                  static_cast<size_t>(position));
          exit(EXIT_FAILURE);
        }
        nexttab[static_cast<int>(cc)]++;
      }
    }
  }
  free(nexttab);
}

template <typename SuftabBaseType,typename SequenceType, bool is_multiseq>
static void gttl_suftab_lightweightcheck(const SequenceType *sequence,
                                         size_t totallength,
                                         [[maybe_unused]] uint8_t wildcard,
                                         [[maybe_unused]] uint8_t padding_char,
                                         const SuftabBaseType *suftab,
                                         size_t numofchars,
                                         [[maybe_unused]]
                                         const size_t *charcount)
{
  size_t count_bits_set = 0;
  size_t previouspos = 0;
  size_t firstspecial = totallength; //NOLINT(misc-const-correctness)
  size_t rangestart = 0;
  uint8_t previouscc = 0;
  Multibitvector<false> startposoccurs(totallength + 1);
  std::vector<Rangewithchar> rangestore;

  rangestore.reserve(numofchars);
  for (size_t idx = 0; idx < totallength; idx++)
  {
    const SuftabBaseType position = suftab[idx];

    if (startposoccurs[static_cast<size_t>(position)])
    {
      fprintf(stderr, "ERROR: suffix with startpos %zu already occurs\n",
              static_cast<size_t>(suftab[idx]));
      exit(EXIT_FAILURE);
    }
    startposoccurs.set(static_cast<size_t>(position));
    count_bits_set++;
    assert(position < totallength);
    const uint8_t cc = static_cast<uint8_t>(sequence[position]);

    if constexpr (is_multiseq)
    {
      if (idx > 0)
      {
        if (isspecial<is_multiseq>(padding_char, cc, wildcard))
        {
          if (firstspecial == totallength)
          {
            firstspecial = idx;
            rangestore.emplace_back(rangestart, idx - 1, previouscc);
          }
          if (isspecial<is_multiseq>(padding_char, previouscc, wildcard))
          {
            if (previouspos > position)
            {
              fprintf(stderr,
                      "incorrect order: %zu = %zu =SPECIAL > "
                      "SPECIAL= %zu  = %zu \n",
                      idx - 1, static_cast<size_t>(position), previouspos, idx);
              exit(EXIT_FAILURE);
            }
          }
        } else
        {
          if (isspecial<is_multiseq>(padding_char, previouscc, wildcard))
          {
            fprintf(stderr,
                    "incorrect order: %zu = %zu =SPECIAL > "
                    "%d= %zu = %zu \n",
                    idx - 1, static_cast<size_t>(position),
                    static_cast<int>(cc),
                    previouspos, idx);
            exit(EXIT_FAILURE);
          } else
          {
            if (previouscc > cc)
            {
              fprintf(stderr, "incorrect order: %zu = %zu=%d > %d=%zu=%zu\n",
                      idx - 1,
                      static_cast<size_t>(position),
                      static_cast<int>(previouscc),
                      static_cast<int>(cc), previouspos, idx);
              exit(EXIT_FAILURE);
            } else
            {
              if (previouscc < cc)
              {
                rangestore.emplace_back(rangestart, idx - 1, previouscc);
                rangestart = idx;
              }
            }
          }
        }
      } else
      {
        if (isspecial<is_multiseq>(padding_char, cc, wildcard))
        {
          firstspecial = 0;
        }
      }
    } else
    {
      if (idx > 0)
      {
        if (previouscc > cc)
        {
          fprintf(stderr, "incorrect order: %zu = %zu=%d > %d=%zu=%zu\n",
                  idx - 1,
                  static_cast<size_t>(position),
                  static_cast<int>(previouscc),
                  static_cast<int>(cc),
                  previouspos,
                  idx);
          exit(EXIT_FAILURE);
        }
        if (previouscc < cc)
        {
          rangestore.emplace_back(rangestart, idx - 1, previouscc);
          rangestart = idx;
        }
      }
    }
    previouscc = cc;
    previouspos = position;
  }
  if (count_bits_set != totallength)
  {
    fprintf(stderr, "ERROR: only %zu of %zu suffixes occur\n", count_bits_set,
            totallength);
    exit(EXIT_FAILURE);
  }
  if (firstspecial == totallength)
  {
    assert(firstspecial > 0);
    rangestore.emplace_back(rangestart, firstspecial - 1, previouscc);
  }
#ifndef NDEBUG
  size_t rangeidx = 0;
  for (size_t charidx = 0; charidx < numofchars; charidx++)
  {
    const size_t count = charcount[charidx];
    if (count > 0)
    {
      assert(rangestore[rangeidx].firstchar_get()
             == static_cast<uint8_t>(charidx));
      assert(rangestore[rangeidx].width_get() == count);
      rangeidx++;
    }
  }
#endif
  gttl_suftab_bk_suffixorder<SuftabBaseType, SequenceType, is_multiseq>
                            (sequence, totallength, wildcard, padding_char,
                             numofchars, suftab, rangestore);
}
#endif
