/*
  Copyright (c) 2013-2022 Stefan Kurtz <stefan.kurtz@uni-hamburg.de>
  Copyright (c) 2013-2022 Center for Bioinformatics, University of Hamburg

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

#ifndef KASAI_HPP
#define KASAI_HPP

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>
#include "utilities/memory_tracker.hpp"
#include "utilities/read_vector.hpp"
#include "inverse_suftab_iter.hpp"

template <typename SuftabBaseType>
static inline std::vector<uint32_t> gttl_lcp13_kasai(
                                      GttlMemoryTracker *memory_tracker,
                                      const uint8_t *sequence,
                                      size_t totallength,
                                      const SuftabBaseType *suftab,
                                      size_t alphasize)
{
  SuftabBaseType *const inverse_suftab = new SuftabBaseType[totallength + 1];

  memory_tracker->track(inverse_suftab,__FILE__,__LINE__,
                        (totallength + 1) * sizeof *inverse_suftab);
  for (size_t idx = 0; idx <= totallength; idx++)
  {
    inverse_suftab[suftab[idx]] = static_cast<SuftabBaseType>(idx);
  }
  std::vector<uint32_t> lcptab;
  lcptab.resize(totallength+1);
  const size_t lcptab_size_in_bytes = (totallength + 1) * sizeof (uint32_t);
  memory_tracker->track(lcptab.data(),__FILE__,__LINE__,lcptab_size_in_bytes);
  lcptab[0] = 0;
  lcptab[totallength] = 0;
  uint32_t lcpvalue = 0;
  for (size_t pos = 0; pos <= totallength; pos++)
  {
    const SuftabBaseType fillpos = inverse_suftab[pos];
    if (fillpos > 0 && static_cast<size_t>(fillpos) < totallength)
    {
      const SuftabBaseType previousstart = suftab[fillpos - 1];
      while (pos + lcpvalue < totallength &&
             previousstart + lcpvalue < totallength)
      {
        const uint8_t cc0 = sequence[pos + lcpvalue];
        const uint8_t cc1 = sequence[previousstart + lcpvalue];
        if (cc0 != cc1 or cc0 >= static_cast<uint8_t>(alphasize))
        {
          break;
        }
        lcpvalue++;
      }
      lcptab[fillpos] = lcpvalue;
    }
    lcpvalue -= (lcpvalue > 0);
  }
  memory_tracker->untrack(inverse_suftab,__FILE__,__LINE__);
  delete[] inverse_suftab;
  return lcptab;
}

template <typename SuftabBaseType>
static inline std::vector<uint32_t> gttl_lcp9_kasai(
                                      GttlMemoryTracker *memory_tracker,
                                      const uint8_t *sequence,
                                      size_t totallength,
                                      const std::string &indexname,
                                      size_t alphasize)
{
  const InverseSuftabReader<SuftabBaseType> inverse_suftab_reader(
                               memory_tracker, indexname, totallength);
  auto inverse_suftab_iter = inverse_suftab_reader.begin();
  std::vector<uint32_t> lcptab;
  lcptab.resize(totallength+1);
  const size_t lcptab_size_in_bytes = (totallength + 1) * sizeof (uint32_t);
  memory_tracker->track(lcptab.data(),__FILE__,__LINE__,lcptab_size_in_bytes);
  lcptab[0] = 0;
  lcptab[totallength] = 0;
  uint32_t lcpvalue = 0;
  const std::string suftab_file{indexname + ".suf"};
  std::vector<SuftabBaseType> suftab
    = gttl_read_vector<SuftabBaseType>(suftab_file);
  memory_tracker->track(suftab.data(),__FILE__,__LINE__,
                        (totallength + 1) * sizeof(SuftabBaseType));
  for (size_t pos = 0; pos <= totallength; pos++)
  {
    const SuftabBaseType fillpos = *inverse_suftab_iter;
    ++inverse_suftab_iter;
    if (fillpos > 0 && static_cast<size_t>(fillpos) < totallength)
    {
      const SuftabBaseType previousstart = suftab[fillpos - 1];
      while (pos + lcpvalue < totallength &&
             previousstart + lcpvalue < totallength)
      {
        const uint8_t cc0 = sequence[pos + lcpvalue];
        const uint8_t cc1 = sequence[previousstart + lcpvalue];
        if (cc0 != cc1 or cc0 >= static_cast<uint8_t>(alphasize))
        {
          break;
        }
        lcpvalue++;
      }
      lcptab[fillpos] = lcpvalue;
    }
    lcpvalue -= (lcpvalue > 0);
  }
  memory_tracker->untrack(suftab.data(),__FILE__,__LINE__);
  return lcptab;
}
#endif
