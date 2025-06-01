/*
  Copyright (c) 2013-2025 Stefan Kurtz <stefan.kurtz@uni-hamburg.de>
  Copyright (c) 2013-2025 Center for Bioinformatics, University of Hamburg

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

#ifndef SK_SAIN_HPP
#define SK_SAIN_HPP

#include <cstdlib>
#include <cstdio>
#include <climits>
#include <cstdint>
#include <iostream>
#include <array>
#include <vector>
#include <map>
#include <type_traits>

#include "utilities/runtime_class.hpp"
#include "utilities/memory_tracker.hpp"
#include "sequences/char_range.hpp"
#include "sequences/char_finder.hpp"
#include "sequences/gttl_multiseq.hpp"
#include "sufcheck.hpp"

#ifdef ACCESS_STATISTICS
static std::map<int,size_t> accesses_per_line;
static std::map<int,size_t> accesses_special_per_line;
static size_t sequential_count = 0;
#endif

class SainbufferKeyValues
{
  public:
  int log_bufsize;
  size_t buf_size,
         cache_size;
  private:
  size_t size;
  bool own_memory;
  public:
  SainbufferKeyValues(size_t sizeof_SuftabBaseType,
                      size_t numofchars,
                      size_t totallength)
    : log_bufsize(std::max(0,21 - (sizeof_SuftabBaseType == size_t(4) ? 1 : 2)
                                - gttl_required_bits(numofchars)))
    , buf_size(size_t(1) << log_bufsize)
    , cache_size(numofchars << log_bufsize)
    , size(sizeof_SuftabBaseType * cache_size +
           sizeof(uint32_t) * numofchars)
    , own_memory(numofchars <= UCHAR_MAX+1 && size * size_t(10) < totallength)
  {
    assert(buf_size <= UINT32_MAX);
  }
  bool has_own_memory(void) const noexcept { return own_memory; }
  size_t size_in_bytes(void) const noexcept
  {
    return size;
  }
};

template<typename SuftabBaseType>
class Sainbuffer
{
  private:
  GttlMemoryTracker *memory_tracker;
  SuftabBaseType *fillptr,
                 *suftab;
  const size_t numofchars;
  const SainbufferKeyValues &sain_buffer_key_values;
  SuftabBaseType *values;
  uint32_t *nextidx;

  public:
  Sainbuffer(GttlMemoryTracker *_memory_tracker,
             const SainbufferKeyValues &_sain_buffer_key_values,
             SuftabBaseType *_suftab,
             SuftabBaseType *_fillptr,
             size_t _numofchars)
    : memory_tracker(_memory_tracker)
    , fillptr(_fillptr)
    , suftab(_suftab)
    , numofchars(_numofchars)
    , sain_buffer_key_values(_sain_buffer_key_values)
    , values(GTTL_TRACK_MALLOC(SuftabBaseType,
                               sain_buffer_key_values.cache_size *
                               sizeof *values))
    , nextidx(GTTL_TRACK_CALLOC(uint32_t,numofchars, sizeof *nextidx))
  { }
  size_t size_in_bytes(void) const noexcept
  {
    return sain_buffer_key_values.size_in_bytes();
  }

  void update(size_t charidx, SuftabBaseType value)
  {
    const size_t offset = charidx << sain_buffer_key_values.log_bufsize;

    values[offset + static_cast<size_t>(nextidx[charidx])] = value;
    if (static_cast<size_t>(nextidx[charidx]) <
        sain_buffer_key_values.buf_size - 1)
    {
      nextidx[charidx]++;
    } else
    {
      const size_t fillptr_offset = fillptr[charidx] - 1;
      SuftabBaseType *writeptr = suftab + fillptr_offset,
                     *readptr = values + offset;
      const SuftabBaseType *readend = readptr + sain_buffer_key_values.buf_size;

      while (readptr < readend)
      {
        *(writeptr--) = *(readptr++);
      }
      nextidx[charidx] = 0;
      fillptr[charidx] -= sain_buffer_key_values.buf_size;
    }
  }

  void flushall(void)
  {
    assert(sain_buffer_key_values.has_own_memory());
    for (size_t charidx = 0; charidx < numofchars; charidx++)
    {
      const SuftabBaseType bufleft
        = static_cast<SuftabBaseType>(nextidx[charidx]);

      if (bufleft > 0)
      {
        SuftabBaseType *writeptr = suftab + fillptr[charidx] - 1,
                       *readptr = values + (charidx <<
                                            sain_buffer_key_values.log_bufsize);
        const SuftabBaseType *readend = readptr + bufleft;

        while (readptr < readend)
        {
          *(writeptr--) = *(readptr++);
        }
        nextidx[charidx] = 0;
        fillptr[charidx] -= bufleft;
      }
    }
  }

  ~Sainbuffer(void)
  {
    GTTL_UNTRACK_ALLOC(values);
    free(values);
    GTTL_UNTRACK_ALLOC(nextidx);
    free(nextidx);
  }
};

class RunTimeAtLevel
{
  std::vector<size_t> run_time_vec;
  public:
  RunTimeAtLevel(void)
  { }
  void add_level(void)
  {
    run_time_vec.push_back(0);
  }
  void set(size_t level, size_t time_microseconds)
  {
    assert(level < run_time_vec.size());
    run_time_vec[level] = time_microseconds;
  }
  void show(FILE *out_fp, size_t level) const
  {
    assert(level < run_time_vec.size());
    size_t sum_upper_level = 0;
    for (size_t idx = level + 1; idx < run_time_vec.size(); idx++)
    {
      sum_upper_level += run_time_vec[idx];
    }
    size_t time_microseconds;
    if (run_time_vec[level] < sum_upper_level)
    {
      time_microseconds = 0;
    } else
    {
      time_microseconds = run_time_vec[level] - sum_upper_level;
    }
    fprintf(out_fp, "# TIME\tlevel %zu (ms):\t%zu\n", level,
            time_microseconds/size_t(1000));
  }
};

static constexpr const char_finder::EncodedNucleotideFinder
                       encoded_nucleotide_finder;
static constexpr const char_finder::EncodeAminoAcidFinder encoded_aa_finder;

typedef enum
{
  GTTL_SAIN_PLAINSEQ,
  GTTL_SAIN_LONGSEQ,
  GTTL_SAIN_MULTISEQ
} GttlSainseqtype;

/* For LONGSEQ, the alphabet is not constant and so we use numofchars, but
   not T_alphasize and so set the latter to 0 */

template <typename SuftabBaseType>
class SuftabResources
{
  GttlMemoryTracker *memory_tracker;
  size_t firstusable;
  const size_t suftabentries;
  SuftabBaseType *suftab;
  public:
  size_t size_in_bytes(void) const noexcept
  {
    return suftabentries * sizeof *suftab;
  }
  SuftabResources(GttlMemoryTracker *_memory_tracker,
                  size_t _suftabentries)
    : memory_tracker(_memory_tracker)
    , firstusable(0)
    , suftabentries(_suftabentries)
    , suftab(GTTL_TRACK_CALLOC(SuftabBaseType,_suftabentries,sizeof *suftab))
  {
    assert (memory_tracker != nullptr);
  }
  void firstusable_set(size_t value)
  {
    assert(firstusable == 0);
    firstusable = value;
  }
  void place_shortest_suffix(void)
  {
    suftab[suftabentries - 1] = suftabentries - 1;
  }
  SuftabBaseType *suftab_ptr_get(void) const
  {
    return suftab;
  }
  auto subtable_assigner(bool do_use_fast_method,size_t number_of_names)
    const noexcept
  {
    std::vector<std::pair<SuftabBaseType *,size_t>> vec_subtable;
    size_t factor = 1;
    for (int idx = 0; idx < (do_use_fast_method ? 3 : 2); idx++)
    {
      const size_t needed = factor * number_of_names;
      if (firstusable + needed <= suftabentries)
      {
        vec_subtable.push_back(std::make_pair(suftab + suftabentries - needed,
                                              0));
      } else
      {
        const size_t this_size = needed * sizeof(SuftabBaseType);
        assert(needed > 0);
        SuftabBaseType *ptr = GTTL_TRACK_MALLOC(SuftabBaseType,this_size);
        vec_subtable.push_back(std::make_pair(ptr,this_size));
        printf("# line %d: extra space[%d] %zu\n",__LINE__,idx,this_size);
      }
      factor *= 2;
    }
    return vec_subtable;
  }
};

template <typename SuftabBaseType,
          GttlSainseqtype T_seqtype, size_t T_alphasize = 0>
class GttlSainseq
{
  static_assert(sizeof(SuftabBaseType) == 4 || sizeof(SuftabBaseType) == 8);
  using Sint = typename std::conditional<sizeof(SuftabBaseType) == 4,
                                         int32_t,
                                         int64_t>::type;
  using T_seq = typename std::conditional
                         <T_seqtype == GTTL_SAIN_PLAINSEQ,
                          uint8_t,
                          typename std::conditional
                                   <T_seqtype == GTTL_SAIN_MULTISEQ,
                                    GttlMultiseq,
                                    SuftabBaseType>::type>::type;
  size_t pos2unique_int(size_t pos) const noexcept
  {
    return pos + UCHAR_MAX + 1;
  }

#define GTTL_SAINUPDATEBUCKETPTR(BUCKETPTR,LASTUPDATECC,CURRENTCC)     \
        if ((BUCKETPTR) != nullptr)                                    \
        {                                                              \
          if ((CURRENTCC) != LASTUPDATECC)                             \
          {                                                            \
            fillptr[LASTUPDATECC] = static_cast<SuftabBaseType>        \
                                               ((BUCKETPTR) - suftab); \
            BUCKETPTR = suftab + fillptr[LASTUPDATECC = CURRENTCC];    \
          }                                                            \
        } else                                                         \
        {                                                              \
          BUCKETPTR = suftab + fillptr[LASTUPDATECC = CURRENTCC];      \
        }

#ifdef SAINSHOWSTATE
  size_t showisStype(const bool *isStype)
  {
    constexpr const bool special_as_pos = T_seqtype == GTTL_SAIN_MULTISEQ;
    for (size_t position = 0; position < totallength; position++)
    {
      std::cout << sainseq_getchar<special_as_pos,true>(position) << " ";
    }
    std::cout << std::endl;
    for (size_t position = 0; position < totallength; position++)
    {
      std::cout << (isStype[position] ? "S" : "L");
    }
    std::cout << "\n ";
    size_t countLMS = 0;
    for (size_t position = 1; position < totallength; position++)
    {
      bool isLMS = (isStype[position] && !isStype[position - 1]) ? true : false;
      std::cout << (isLMS ? "*" : " ");
      if (isLMS)
      {
        countLMS++;
      }
    }
    std::cout << std::endl;
    return countLMS;
  }
#endif

  private:
  GttlMemoryTracker *memory_tracker;
  const size_t totallength,
               numofchars,
               nonspecialentries;
#define WITH_BRANCHING
#ifndef WITH_BRANCHING
#define WITH_BRANCHINGconst /* Nothing */
  std::array<size_t,T_alphasize + 1> return_values;
#else
#define WITH_BRANCHINGconst const
#endif
  const size_t *charcount_arr;
  SuftabBaseType currentround,
                 *bucketsizeptr,
                 *bucketfillptr,
                 *roundtable, /* only for multiseq/plain */
                 *sstarfirstcharcount;
  const T_seq *seq;
  bool bucketsizepoints2suftab,
       bucketfillptrpoints2suftab,
       roundtablepoints2suftab;
  const char *multiseq_seq_ptr;
  size_t this_size;

  template <bool special_as_pos,bool sequential = false>
  size_t sainseq_getchar(size_t position,GTTL_UNUSED int line)
    WITH_BRANCHINGconst
    noexcept
  {
    assert(position < totallength);
#ifdef ACCESS_STATISTICS
    accesses_per_line[line]++;
    sequential_count += sequential;
#endif
    if constexpr (T_seqtype == GTTL_SAIN_MULTISEQ)
    {
#ifdef WITH_BRANCHING
      const size_t cc = static_cast<size_t>(multiseq_seq_ptr[position]);
      if constexpr (special_as_pos)
      {
#ifdef ACCESS_STATISTICS
        accesses_special_per_line[line] += (cc == T_alphasize);
#endif
        return cc < T_alphasize ? cc : pos2unique_int(position);
      } else
      {
        return cc;
      }
#else
      if constexpr (special_as_pos)
      {
        return_values[T_alphasize] = pos2unique_int(position);
      }
      const size_t cc = static_cast<size_t>(multiseq_seq_ptr[position]);
      return return_values[cc];
#endif
    } else
    {
      return seq[position];
    }
    return 0;
  }

  bool use_fast_method(void) const noexcept
  {
    constexpr const int word_size = sizeof(SuftabBaseType) * CHAR_BIT;
    constexpr const SuftabBaseType first2bits
      = static_cast<SuftabBaseType>(3) << (word_size - 2); /* \(11^{w-2}\) */
    const size_t maxvalue = totallength + 1;
    return maxvalue < first2bits &&
           totallength > size_t(1024) &&
           totallength >= 2 * numofchars;
  }

  void sain_endbuckets(void)
  {
    SuftabBaseType previous = bucketfillptr[0] = bucketsizeptr[0];
    for (size_t charidx = size_t(1); charidx < numofchars; charidx++)
    {
      previous += bucketsizeptr[charidx];
      bucketfillptr[charidx] = previous;
    }
  }

  size_t sain_insertSstarsuffixes(SuftabBaseType *suftab)
  {
    size_t nextcc = pos2unique_int(totallength), countSstartype = 0;
    SuftabBaseType *fillptr = bucketfillptr;
    bool nextisStype = true;
#ifdef SAINSHOWSTATE
    size_t countLMS;
    bool *isStype = GTTL_TRACK_MALLOC(bool,totallength * sizeof *isStype);
#endif

    sain_endbuckets();
    for (size_t position = totallength - 1; /* Nothing */; position--)
    {
      constexpr const bool special_as_pos = T_seqtype == GTTL_SAIN_MULTISEQ;
      const size_t currentcc = sainseq_getchar<special_as_pos,true>
                                              (position,__LINE__);
      const bool currentisStype
        = currentcc < nextcc || (currentcc == nextcc && nextisStype);
#ifdef SAINSHOWSTATE
      isStype[position] = currentisStype;
#endif
      if (!currentisStype && nextisStype)
      {
        countSstartype++;
        if constexpr (T_seqtype != GTTL_SAIN_LONGSEQ)
        {
          sstarfirstcharcount[nextcc]++;
        }
        suftab[--fillptr[nextcc]] = position;
#ifdef SAINSHOWSTATE
        printf("Sstar.suftab[%zu]=%zu\n", static_cast<size_t>(fillptr[nextcc]),
               position + 1);
#endif
      }
      nextisStype = currentisStype;
      nextcc = currentcc;

      if (position == 0)
      {
        break;
      }
    }
#ifdef SAINSHOWSTATE
    countLMS = showisStype(isStype);
    printf("countLMS=%zu\n", countLMS);
    free(isStype);
#endif
    assert(2 * countSstartype <= totallength);
    return countSstartype;
  }

  size_t sain_buffered_insertSstarsuffixes(const SainbufferKeyValues
                                             &sain_buffer_key_values,
                                           SuftabBaseType *suftab)
  {
    SuftabBaseType *fillptr = bucketfillptr;

    Sainbuffer<SuftabBaseType> sainbuffer
                               (memory_tracker,
                                sain_buffer_key_values,
                                suftab, fillptr, numofchars);
    bool nextisStype = true;
#ifdef SAINSHOWSTATE
    size_t countLMS;
    bool *isStype = GTTL_TRACK_MALLOC(bool,totallength * sizeof *isStype);
#endif

    sain_endbuckets();
    size_t nextcc = pos2unique_int(totallength), countSstartype = 0;
    for (size_t position = totallength - 1; /* Nothing */; position--)
    {
      constexpr const bool special_as_pos = T_seqtype == GTTL_SAIN_MULTISEQ;
      const size_t currentcc = sainseq_getchar<special_as_pos,true>
                                              (position,__LINE__);
      const bool currentisStype
        = currentcc < nextcc || (currentcc == nextcc && nextisStype);
#ifdef SAINSHOWSTATE
      isStype[position] = currentisStype;
#endif
      if (!currentisStype && nextisStype)
      {
        countSstartype++;
        if constexpr (T_seqtype != GTTL_SAIN_LONGSEQ)
        {
          sstarfirstcharcount[nextcc]++;
        }
        sainbuffer.update(nextcc, position);

#ifdef SAINSHOWSTATE
        printf("Sstar.suftab[%zu]=%zu\n", static_cast<size_t>(fillptr[nextcc]),
               position + 1);
#endif
      }
      nextisStype = currentisStype;
      nextcc = currentcc;

      if (position == 0)
      {
        break;
      }
    }
#ifdef SAINSHOWSTATE
    countLMS = showisStype(isStype);
    printf("countLMS=%zu\n",countLMS);
    free(isStype);
#endif
    sainbuffer.flushall();
    assert(2 * countSstartype <= totallength);
    return countSstartype;
  }

  void sain_incrementfirstSstar(SuftabBaseType *suftab)
  {
    SuftabBaseType sum = 0;

    assert(roundtable != nullptr);
    for (size_t charidx = 0; charidx < numofchars; charidx++)
    {
      sum += bucketsizeptr[charidx];
      assert(bucketfillptr[charidx] <= sum);
      if (bucketfillptr[charidx] < sum)
      {
        suftab[bucketfillptr[charidx]] += totallength;
      }
      roundtable[charidx] = 0;
      roundtable[charidx + numofchars] = 0;
    }
  }

  void sain_startbuckets(void)
  {
    SuftabBaseType previous = bucketfillptr[0] = 0;
    for (size_t charidx = size_t(1); charidx < numofchars; charidx++)
    {
      previous += bucketsizeptr[charidx - 1];
      bucketfillptr[charidx] = previous;
    }
  }

  void sain_PLAIN_MULTISEQ_induceLtypesuffixes1(Sint *suftab)
  {
    assert(numofchars == T_alphasize);
    size_t lastupdatecc = 0;
    SuftabBaseType *fillptr = bucketfillptr;
    Sint *bucketptr = nullptr;

    for (size_t suftab_idx = 0; suftab_idx < nonspecialentries; suftab_idx++)
    {
      Sint position = suftab[suftab_idx];
      if (position > 0)
      {
        constexpr const bool special_as_pos = T_seqtype == GTTL_SAIN_MULTISEQ;
        const size_t currentcc = sainseq_getchar<special_as_pos>
                                                (position,__LINE__);
        assert(T_seqtype == GTTL_SAIN_MULTISEQ || currentcc < T_alphasize);
        if (currentcc < T_alphasize)
        {
          GTTL_SAINUPDATEBUCKETPTR(bucketptr,lastupdatecc,currentcc);
          /* negative => position does not derive L-suffix
          positive => position may derive L-suffix */
          assert(suftab + suftab_idx < bucketptr);
          position--;
          const size_t leftcontextcc = sainseq_getchar<special_as_pos>
                                                      (position,__LINE__);
          *bucketptr++ = (leftcontextcc < currentcc) ? ~position : position;
          suftab[suftab_idx] = 0;
#ifdef SAINSHOWSTATE
          printf("line %d: L-induce: suftab[%zu]=%d\n", __LINE__,
                 static_cast<size_t>(bucketptr - 1 - suftab),
                 *(bucketptr - 1));
#endif
        } else
        {
          suftab[suftab_idx] = 0;
        }
      } else
      {
        if (position < 0)
        {
          suftab[suftab_idx] = ~position;
        }
      }
    }
  }

  void sain_PLAIN_MULTISEQ_fast_induceLtypesuffixes1(Sint *suftab)
  {
    assert(numofchars == T_alphasize);
    size_t lastupdatecc = 0;
    SuftabBaseType *fillptr = bucketfillptr;
    Sint *bucketptr = nullptr;

    assert(roundtable != nullptr);
    currentround = 0;
    for (size_t suftab_idx = 0; suftab_idx < nonspecialentries; suftab_idx++)
    {
      Sint position = suftab[suftab_idx];
      if (position > 0)
      {
        if (position >= static_cast<Sint>(totallength))
        {
          currentround++;
          position -= static_cast<Sint>(totallength);
        }
        /* false as used for PLAINSEQ */
        constexpr const bool special_as_pos = T_seqtype == GTTL_SAIN_MULTISEQ;
        const size_t currentcc = sainseq_getchar<special_as_pos>
                                                (position,__LINE__);
        assert(T_seqtype == GTTL_SAIN_MULTISEQ || currentcc < T_alphasize);
        if (currentcc < T_alphasize)
        {
          if (position > 0)
          {
            position--;
            const size_t leftcontextcc = sainseq_getchar<special_as_pos>
                                                        (position,__LINE__);
            const size_t t = (currentcc << 1) |
                             (leftcontextcc < currentcc ? size_t(1) : 0);
            assert(roundtable[t] <= currentround);
            if (roundtable[t] < currentround)
            {
              position += static_cast<long>(totallength);
              roundtable[t] = currentround;
            }
            GTTL_SAINUPDATEBUCKETPTR(bucketptr,lastupdatecc,currentcc);
            /* negative => position does not derive L-suffix
               positive => position may derive L-suffix */
            assert(suftab + suftab_idx < bucketptr);
            *bucketptr++ = (t & size_t(1)) ? ~position : position;
            suftab[suftab_idx] = 0;
#ifdef SAINSHOWSTATE
            printf("line %d: L-induce: suftab[%zu]=%d\n", __LINE__,
                   static_cast<size_t>(bucketptr - 1 - suftab),
                   *(bucketptr - 1));
#endif
          }
        } else
        {
          suftab[suftab_idx] = 0;
        }
      } else
      {
        if (position < 0)
        {
          suftab[suftab_idx] = ~position;
        }
      }
    }
  }

  void sain_LONGSEQ_induceLtypesuffixes1(Sint *suftab)
  {
    size_t lastupdatecc = 0;
    SuftabBaseType *fillptr = bucketfillptr;
    Sint *bucketptr = nullptr;

    assert(roundtable == nullptr);
    for (size_t suftab_idx = 0; suftab_idx < nonspecialentries; suftab_idx++)
    {
      Sint position = suftab[suftab_idx];
      if (position > 0)
      {
        /* false as used for LONGSEQ */
        const size_t currentcc = sainseq_getchar<false>(position,__LINE__);
        position--;
        const size_t leftcontextcc = sainseq_getchar<false>(position,__LINE__);
        GTTL_SAINUPDATEBUCKETPTR(bucketptr,lastupdatecc,currentcc);
        /* negative => position does not derive L-suffix
          positive => position may derive L-suffix */
        assert(suftab + suftab_idx < bucketptr);
        *bucketptr++ = (leftcontextcc < currentcc) ? ~position : position;
        suftab[suftab_idx] = 0;
#ifdef SAINSHOWSTATE
        printf("line %d: L-induce: suftab[%zu]=%d\n", __LINE__,
               static_cast<size_t>(bucketptr - 1 - suftab),
               *(bucketptr - 1));
#endif
      } else
      {
        if (position < 0)
        {
          suftab[suftab_idx] = ~position;
        }
      }
    }
  }

  void sain_LONGSEQ_fast_induceLtypesuffixes1(Sint *suftab)
  {
    size_t lastupdatecc = 0;
    SuftabBaseType *fillptr = bucketfillptr;
    Sint *bucketptr = nullptr;

    assert(roundtable != nullptr);
    currentround = 0;
    for (size_t suftab_idx = 0; suftab_idx < nonspecialentries; suftab_idx++)
    {
      Sint position = suftab[suftab_idx];
      if (position > 0)
      {
        if (position >= static_cast<Sint>(totallength))
        {
          currentround++;
          position -= static_cast<Sint>(totallength);
        }
        /* false as used for LONGSEQ */
        const size_t currentcc = sainseq_getchar<false>(position,__LINE__);
        if (position > 0)
        {
          position--;
          const size_t leftcontextcc = sainseq_getchar<false>(position,
                                                              __LINE__);
          const size_t t = (currentcc << 1) |
                           (leftcontextcc < currentcc ? size_t(1) : 0);
          assert(currentcc > 0 && roundtable[t] <= currentround);
          if (roundtable[t] < currentround)
          {
            position += static_cast<Sint>(totallength);
            roundtable[t] = currentround;
          }
          GTTL_SAINUPDATEBUCKETPTR(bucketptr,lastupdatecc,currentcc);
          /* negative => position does not derive L-suffix
            positive => position may derive L-suffix */
          assert(suftab + suftab_idx < bucketptr);
          *bucketptr++ = (t & size_t(1)) ? ~position : position;
          suftab[suftab_idx] = 0;
#ifdef SAINSHOWSTATE
          printf("line %d: L-induce: suftab[%zu]=%d\n", __LINE__,
                 static_cast<size_t>(bucketptr - 1 - suftab), *(bucketptr - 1));
#endif
        }
      } else
      {
        if (position < 0)
        {
          suftab[suftab_idx] = ~position;
        }
      }
    }
  }

  void sain_induceLtypesuffixes1(Sint *suftab)
  {
    if constexpr (T_seqtype == GTTL_SAIN_LONGSEQ)
    {
      if (roundtable == nullptr)
      {
        sain_LONGSEQ_induceLtypesuffixes1(suftab);
      } else
      {
        sain_LONGSEQ_fast_induceLtypesuffixes1(suftab);
      }
    } else
    {
      if (roundtable == nullptr)
      {
        sain_PLAIN_MULTISEQ_induceLtypesuffixes1(suftab);
      } else
      {
        sain_PLAIN_MULTISEQ_fast_induceLtypesuffixes1(suftab);
      }
    }
  }

  void sain_special_singleSinduction1(Sint *suftab, Sint position)
  {
    constexpr const bool special_as_pos = T_seqtype == GTTL_SAIN_MULTISEQ;
    /* position is always a position before a range of special symbols
       or the position is totallength - 1. The number is accesses is
       equal to the number of ranges + 1 */
    const size_t currentcc = sainseq_getchar<special_as_pos>
                                            (position,__LINE__);
    if (currentcc < numofchars)
    {
      const size_t putidx = static_cast<size_t>(--bucketfillptr[currentcc]);
      assert(position > 0);
      position--;
      /* The number of these accesses is upper bounded by the number
         above */
      const size_t leftcontextcc = sainseq_getchar<special_as_pos>
                                                  (position,__LINE__);
      if (roundtable != nullptr)
      {
        const size_t t
          = (currentcc << 1) | (leftcontextcc > currentcc ? size_t(1) : 0);
        assert(roundtable[t] <= currentround);
        if (roundtable[t] < currentround)
        {
          roundtable[t] = currentround;
        }
        position += totallength;
      }
      suftab[putidx] = (leftcontextcc > currentcc) ? ~(position + 1) : position;
#ifdef SAINSHOWSTATE
      printf("end S-induce: suftab[%zu]=%ld\n", putidx,
             static_cast<long>(suftab[putidx]));
#endif
    }
  }

  void induceStypes1fromspecialranges(Sint *suftab)
  {
    // get the position of the padding chars and wildcards from second last to
    // first
    static_assert(T_seqtype == GTTL_SAIN_MULTISEQ);
    assert(seq->sequences_number_get() > 0);
    constexpr const bool forward = false, invert = true;

    if constexpr (T_alphasize == size_t(4))
    {
      using ThisCharRange
        = GttlCharRange<char_finder::EncodedNucleotideFinder,
                        encoded_nucleotide_finder, forward, invert>;
      assert(multiseq_seq_ptr != nullptr);
      ThisCharRange ranger(multiseq_seq_ptr, totallength);
      for (auto &&range : ranger)
      {
        if (std::get<0>(range) > 1)
        {
          sain_special_singleSinduction1(suftab,
                                         static_cast<Sint>(std::get<0>(range)
                                                             - 1));
        }
      }
    } else
    {
      static_assert(T_alphasize == size_t(20));
      assert(multiseq_seq_ptr != nullptr);
      using ThisCharRange = GttlCharRange<char_finder::EncodeAminoAcidFinder,
                                          encoded_aa_finder, forward, invert>;
      ThisCharRange ranger(multiseq_seq_ptr, totallength);
      for (auto &&range : ranger)
      {
        if (std::get<0>(range) > 1)
        {
          sain_special_singleSinduction1(suftab,
                                         static_cast<Sint>(std::get<0>(range)
                                                             - 1));
        }
      }
    }
  }

  void sain_PLAIN_MULTISEQ_induceStypesuffixes1(Sint *suftab)
  {
    assert(numofchars == T_alphasize);
    size_t lastupdatecc = 0;
    SuftabBaseType *fillptr = bucketfillptr;
    Sint *bucketptr = nullptr;
    assert(roundtable == nullptr);
    sain_special_singleSinduction1(suftab, static_cast<Sint>(totallength - 1));
    if constexpr (T_seqtype == GTTL_SAIN_MULTISEQ)
    {
      induceStypes1fromspecialranges(suftab);
    }
    for (size_t suftab_idx = nonspecialentries; suftab_idx > 0; /* Nothing */)
    {
      suftab_idx--;
      Sint position = suftab[suftab_idx];
      if (position > 0)
      {
        constexpr const bool special_as_pos = T_seqtype == GTTL_SAIN_MULTISEQ;
        const size_t currentcc = sainseq_getchar<special_as_pos>
                                                (position,__LINE__);
        assert(T_seqtype == GTTL_SAIN_MULTISEQ || currentcc < T_alphasize);
        if (currentcc < T_alphasize)
        {
          GTTL_SAINUPDATEBUCKETPTR(bucketptr,lastupdatecc,currentcc);
          assert(bucketptr - 1 < suftab + suftab_idx);
          position--;
          const size_t leftcontextcc = sainseq_getchar<special_as_pos>
                                                      (position,__LINE__);
          *(--bucketptr) = (leftcontextcc > currentcc) ? ~(position + 1)
                                                       : position;

#ifdef SAINSHOWSTATE
          printf("line %d: S-induce: suftab[%zu]=%d\n", __LINE__,
                 static_cast<size_t>(bucketptr - suftab), *bucketptr);
#endif
        }
        suftab[suftab_idx] = 0;
      }
    }
  }

  void sain_PLAIN_MULTISEQ_fast_induceStypesuffixes1(Sint *suftab)
  {
    assert(numofchars == T_alphasize);
    size_t lastupdatecc = 0;
    SuftabBaseType *fillptr = bucketfillptr;
    Sint *bucketptr = nullptr;

    assert(roundtable != nullptr);
    sain_special_singleSinduction1(suftab, static_cast<Sint>(totallength - 1));
    if constexpr (T_seqtype == GTTL_SAIN_MULTISEQ)
    {
      induceStypes1fromspecialranges(suftab);
    }
    for (size_t suftab_idx = nonspecialentries; suftab_idx > 0; /* Nothing */)
    {
      suftab_idx--;
      Sint position = suftab[suftab_idx];
      if (position > 0)
      {
        /* we use the signed position values and only the positive have to
           be further processed here. So this reduces the number of sequence
           positions to represent by 50%. As we also sometime need to
           increment the positions by totallength this reduces the number
           of sequence position by another 50%. So we only can represent
           25% of the sequence position which could be represented by a
           value of SuftabBaseType. If this is uint32_t, then we can represent
           positions upto 2^{32}/4=2^{30}. */
        if (position >= static_cast<Sint>(totallength))
        {
          currentround++;
          position -= static_cast<Sint>(totallength);
        }
        if (position > 0)
        {
          const size_t currentcc = sainseq_getchar<false>(position,__LINE__);
          /* It seems that currentcc is never a special symbol */
          position--;
          constexpr const bool special_as_pos = T_seqtype == GTTL_SAIN_MULTISEQ;
          const size_t leftcontextcc = sainseq_getchar<special_as_pos>
                                                      (position,__LINE__);
          const size_t t = (currentcc << 1) |
                           (leftcontextcc > currentcc ? size_t(1) : 0);
          assert(roundtable[t] <= currentround);
          if (roundtable[t] < currentround)
          {
            position += totallength;
            roundtable[t] = currentround;
          }
          GTTL_SAINUPDATEBUCKETPTR(bucketptr,lastupdatecc,currentcc);
          assert(bucketptr - 1 < suftab + suftab_idx);
          *(--bucketptr) = (t & size_t(1)) ? ~(position + 1) : position;
#ifdef SAINSHOWSTATE
          printf("line %d: S-induce: suftab[%zu]=%d\n", __LINE__,
                 static_cast<size_t>(bucketptr - suftab), *bucketptr);
#endif
        }
        suftab[suftab_idx] = 0;
      }
    }
  }

  void sain_LONGSEQ_induceStypesuffixes1(Sint *suftab)
  {
    assert(T_alphasize == 0);
    size_t lastupdatecc = 0;
    SuftabBaseType *fillptr = bucketfillptr;
    Sint *bucketptr = nullptr;

    assert(roundtable == nullptr);
    sain_special_singleSinduction1(suftab, static_cast<Sint>(totallength - 1));
    for (size_t suftab_idx = nonspecialentries; suftab_idx > 0; /* Nothing */)
    {
      suftab_idx--;
      Sint position = suftab[suftab_idx];
      if (position > 0)
      {
        /* false as used for LONGSEQ */
        const size_t currentcc = sainseq_getchar<false>(position,__LINE__);
        position--;
        const size_t leftcontextcc = sainseq_getchar<false>(position,__LINE__);
        GTTL_SAINUPDATEBUCKETPTR(bucketptr,lastupdatecc,currentcc);
        assert(bucketptr - 1 < suftab + suftab_idx);
        *(--bucketptr) = (leftcontextcc > currentcc) ? ~(position + 1)
                                                     : position;
#ifdef SAINSHOWSTATE
        printf("line %d: S-induce: suftab[%zu]=%d\n", __LINE__,
               static_cast<size_t>(bucketptr - suftab), *bucketptr);
#endif
        suftab[suftab_idx] = 0;
      }
    }
  }

  void sain_LONGSEQ_fast_induceStypesuffixes1(Sint *suftab)
  {
    assert(T_alphasize == 0);
    size_t lastupdatecc = 0;
    SuftabBaseType *fillptr = bucketfillptr;
    Sint *bucketptr = nullptr;

    assert(roundtable != nullptr);
    sain_special_singleSinduction1(suftab, static_cast<Sint>(totallength - 1));
    for (size_t suftab_idx = nonspecialentries; suftab_idx > 0; /* Nothing */)
    {
      suftab_idx--;
      Sint position = suftab[suftab_idx];
      if (position > 0)
      {
        if (position >= static_cast<Sint>(totallength))
        {
          currentround++;
          position -= static_cast<Sint>(totallength);
        }
        if (position > 0)
        {
          /* false as used for LONGSEQ */
          const size_t currentcc = sainseq_getchar<false>(position,__LINE__);

          position--;
          const size_t leftcontextcc = sainseq_getchar<false>(position,
                                                              __LINE__);
          const size_t t = (currentcc << 1) |
                           (leftcontextcc > currentcc ? size_t(1) : 0);
          assert(roundtable[t] <= currentround);
          if (roundtable[t] < currentround)
          {
            position += totallength;
            roundtable[t] = currentround;
          }
          GTTL_SAINUPDATEBUCKETPTR(bucketptr,lastupdatecc,currentcc);
          assert(bucketptr - 1 < suftab + suftab_idx);
          *(--bucketptr) = (t & size_t(1)) ? ~(position + 1) : position;
#ifdef SAINSHOWSTATE
          printf("line %d, S-induce: suftab[%zu]=%d\n", __LINE__,
                 static_cast<size_t>(bucketptr - suftab), *bucketptr);
#endif
        }
        suftab[suftab_idx] = 0;
      }
    }
  }

  void sain_induceStypesuffixes1(Sint *suftab)
  {
    if constexpr (T_seqtype == GTTL_SAIN_LONGSEQ)
    {
      assert(T_alphasize == 0);
      if (roundtable == nullptr)
      {
        sain_LONGSEQ_induceStypesuffixes1(suftab);
      } else
      {
        sain_LONGSEQ_fast_induceStypesuffixes1(suftab);
      }
    } else
    {
      if (roundtable == nullptr)
      {
        sain_PLAIN_MULTISEQ_induceStypesuffixes1(suftab);
      } else
      {
        sain_PLAIN_MULTISEQ_fast_induceStypesuffixes1(suftab);
      }
    }
  }

  void sain_assignSstarlength(SuftabBaseType *lentab)
  {
    bool nextisStype = true;
    size_t nextSstartypepos = totallength,
           nextcc = pos2unique_int(totallength);

    for (size_t position = totallength - 1; /* Nothing */; position--)
    {
      constexpr const bool special_as_pos = T_seqtype == GTTL_SAIN_MULTISEQ;
      const size_t currentcc = sainseq_getchar<special_as_pos,true>
                                              (position,__LINE__);
      const bool currentisStype
        = currentcc < nextcc || (currentcc == nextcc && nextisStype);
      if (!currentisStype && nextisStype)
      {
        assert(position < nextSstartypepos);
        lentab[(position + 1)/2]
          = static_cast<SuftabBaseType>(nextSstartypepos - position);
        nextSstartypepos = position + 1;
      }
      nextisStype = currentisStype;
      nextcc = currentcc;
      if (position == 0)
      {
        break;
      }
    }
  }

  int sain_compare_Sstarstrings(size_t start1, size_t start2, size_t len)
  {
    const size_t end1 = start1 + len;

    assert(start1 <= totallength && start2 <= totallength && start1 != start2);
    while (start1 < end1)
    {
      size_t cc1, cc2;

      if (start1 == totallength)
      {
        assert(start1 > start2);
        return 1;
      }
      if (start2 == totallength)
      {
        assert(start1 <= start2);
        return -1;
      }
      constexpr const bool special_as_pos = T_seqtype == GTTL_SAIN_MULTISEQ;
      cc1 = sainseq_getchar<special_as_pos>(start1,__LINE__);
      cc2 = sainseq_getchar<special_as_pos>(start2,__LINE__);
      if (cc1 < cc2)
      {
        return -1;
      }
      if (cc1 > cc2)
      {
        return 1;
      }
      start1++;
      start2++;
    }
    return 0;
  }

  size_t sain_assignSstarnames(size_t countSstartype, SuftabBaseType *suftab)
  {
    SuftabBaseType *secondhalf = suftab + countSstartype,
                   previouspos = suftab[0],
                   previouslen = secondhalf[previouspos/2];
    size_t currentname = static_cast<size_t>(1);

    secondhalf[previouspos/2] = static_cast<SuftabBaseType>(currentname);
    for (size_t suftab_idx = size_t(1); suftab_idx < countSstartype;
         suftab_idx++)
    {
      const SuftabBaseType position = suftab[suftab_idx],
                           currentlen = secondhalf[position/2];

      currentname += (previouslen != currentlen) ||
                     (sain_compare_Sstarstrings(previouspos,position,currentlen)
                      == -1);
      /* write the names in order of positions. As the positions of
        the Sstar suffixes differ by at least 2, the used address
        is unique */
      previouslen = currentlen;
      secondhalf[position/2] = static_cast<SuftabBaseType>(currentname);
      previouspos = position;
    }
    return currentname;
  }

  void sain_expandorder2original(size_t numberofsuffixes,
                                 SuftabBaseType *suftab)
  {
    SuftabBaseType *sstarsuffixes = suftab + 2 * numberofsuffixes;
    size_t nextcc = pos2unique_int(totallength);
    bool nextisStype = true;
    if constexpr (T_seqtype == GTTL_SAIN_LONGSEQ)
    {
      assert(T_alphasize == 0 && sstarfirstcharcount == nullptr);
      sstarfirstcharcount = bucketfillptr;
      for (size_t charidx = 0; charidx < numofchars; charidx++)
      {
        sstarfirstcharcount[charidx] = 0;
        bucketsizeptr[charidx] = 0;
      }
    }
    for (SuftabBaseType position = totallength - 1; /* Nothing */; position--)
    {
      constexpr const bool special_as_pos = T_seqtype == GTTL_SAIN_MULTISEQ;
      const size_t currentcc
        = sainseq_getchar<special_as_pos,true>(position,__LINE__);
      const bool currentisStype
        = currentcc < nextcc || (currentcc == nextcc && nextisStype);
      if (!currentisStype && nextisStype)
      {
        if constexpr (T_seqtype == GTTL_SAIN_LONGSEQ)
        {
          sstarfirstcharcount[nextcc]++;
        }
        *(--sstarsuffixes) = position + 1;
      }
      if constexpr (T_seqtype == GTTL_SAIN_LONGSEQ)
      {
        bucketsizeptr[currentcc]++;
      }
      nextisStype = currentisStype;
      nextcc = currentcc;
      if (position == 0)
      {
        break;
      }
    }
    for (size_t suftab_idx = 0; suftab_idx < numberofsuffixes; suftab_idx++)
    {
      if constexpr (T_seqtype != GTTL_SAIN_LONGSEQ)
      {
        assert(suftab[suftab_idx] < numberofsuffixes);
      }
      suftab[suftab_idx] = sstarsuffixes[suftab[suftab_idx]];
    }
  }

  void sain_determineSstarfirstchardist(void)
  {
    size_t nextcc = pos2unique_int(totallength);
    bool nextisStype = true;

    assert(totallength > 0 && T_seqtype == GTTL_SAIN_LONGSEQ);
    for (size_t pos = totallength - 1; /* Nothing */; pos--)
    {
      /* false as use for LONGSEQ */
      const size_t currentcc = sainseq_getchar<false,true>(pos,__LINE__);
      const bool currentisStype
        = currentcc < nextcc || (currentcc == nextcc && nextisStype);
      if (!currentisStype && nextisStype)
      {
        sstarfirstcharcount[nextcc]++;
      }
      bucketsizeptr[currentcc]++;
      nextisStype = currentisStype;
      nextcc = currentcc;
      if (pos == 0)
      {
        break;
      }
    }
  }

  int sain_compare_suffixes(size_t start1, size_t start2)
    WITH_BRANCHINGconst noexcept
  {
    assert(start1 <= totallength && start2 <= totallength && start1 != start2);
    constexpr const bool special_as_pos = T_seqtype == GTTL_SAIN_MULTISEQ;

    while (true)
    {
      if (start1 == totallength)
      {
        assert(start1 > start2);
        return 1;
      }
      if (start2 == totallength)
      {
        assert(start1 < start2);
        return -1;
      }
      const size_t cc1 = sainseq_getchar<special_as_pos>(start1,__LINE__);
      const size_t cc2 = sainseq_getchar<special_as_pos>(start2,__LINE__);
      if (cc1 < cc2)
      {
        return -1;
      }
      if (cc1 > cc2)
      {
        return 1;
      }
      start1++;
      start2++;
    }
  }

  void sain_checkorder(const SuftabBaseType *suftab, size_t start, size_t end)
    WITH_BRANCHINGconst noexcept
  {
    for (size_t idx = start + 1; idx <= end; idx++)
    {
      int cmp = sain_compare_suffixes(suftab[idx - 1], suftab[idx]);

      if (cmp != -1)
      {
        fprintf(stderr,
                "%s: check interval [%zu,%zu] at idx=%zu: suffix "
                "%zu >= %zu\n",
                __func__, start, end, idx, static_cast<size_t>(suftab[idx - 1]),
                static_cast<size_t>(suftab[idx]));
        exit(EXIT_FAILURE);
      }
    }
  }

  void sain_insertsortedSstarsuffixes(SuftabBaseType *suftab, size_t readidx)
  {
    size_t fillidx = nonspecialentries;

    for (size_t cc = numofchars - 1; /* Nothing */; cc--)
    {
      if (sstarfirstcharcount[cc] > 0)
      {
        const size_t putidx = fillidx - 1;

        assert(readidx <= putidx);
        if (readidx < putidx)
        {
          for (SuftabBaseType offset = 0;
               offset < sstarfirstcharcount[cc]; offset++)
          {
            suftab[putidx - offset] = suftab[readidx - offset];
            suftab[readidx - offset] = 0;
#ifdef SAINSHOWSTATE
            printf("insertsorted: suftab[%zu]=%zu\n", putidx - offset,
                   static_cast<size_t>(suftab[putidx - offset]));
            printf("insertsorted: suftab[%u]=undef\n", readidx - offset);
#endif
          }
        }
      }
      assert(fillidx >= bucketsizeptr[cc] &&
             bucketsizeptr[cc] >= sstarfirstcharcount[cc]);
      fillidx -= bucketsizeptr[cc];
      if (bucketsizeptr[cc] > sstarfirstcharcount[cc])
      {
        sain_set_undefined<false>(suftab, fillidx,
                                  fillidx + bucketsizeptr[cc]
                                          - sstarfirstcharcount[cc] - 1);
      }
      readidx -= sstarfirstcharcount[cc];
      if (cc == 0)
      {
        break;
      }
    }
  }

  void sain_PLAIN_MULTISEQ_induceLtypesuffixes2(Sint *suftab)
  {
    assert(T_alphasize == numofchars);
    SuftabBaseType *fillptr = bucketfillptr;
    size_t lastupdatecc = 0;
    Sint *bucketptr = nullptr;

    for (size_t suftab_idx = 0; suftab_idx < nonspecialentries; suftab_idx++)
    {
      Sint position = suftab[suftab_idx];
      suftab[suftab_idx] = ~position;
      if (position > 0)
      {
        position--;
        constexpr const bool special_as_pos = T_seqtype == GTTL_SAIN_MULTISEQ;
        const size_t currentcc = sainseq_getchar<special_as_pos>(position,
                                                                 __LINE__);
        assert(T_seqtype == GTTL_SAIN_MULTISEQ || currentcc < T_alphasize);
        if (currentcc < T_alphasize)
        {
          assert(currentcc > 0);
          GTTL_SAINUPDATEBUCKETPTR(bucketptr,lastupdatecc,currentcc)
          assert(suftab + suftab_idx < bucketptr);
          *bucketptr++ = (position > 0 &&
                          sainseq_getchar<special_as_pos>(position - 1,__LINE__)
                            < lastupdatecc)
                             ? ~position
                             : position;

#ifdef SAINSHOWSTATE
          assert(bucketptr != nullptr);
          printf("line %d: L-induce: suftab[%zu]=%d\n", __LINE__,
                 static_cast<size_t>(bucketptr - 1 - suftab), *(bucketptr - 1));
#endif
        }
      }
    }
  }

  void sain_LONGSEQ_induceLtypesuffixes2(Sint *suftab)
  {
    assert(T_alphasize == 0);
    size_t lastupdatecc = 0;
    SuftabBaseType *fillptr = bucketfillptr;
    Sint *bucketptr = nullptr;
    for (size_t suftab_idx = 0; suftab_idx < nonspecialentries; suftab_idx++)
    {
      Sint position = suftab[suftab_idx];
      suftab[suftab_idx] = ~position;
      if (position > 0)
      {
        position--;
        /* false as used for LONGSEQ */
        const size_t currentcc = sainseq_getchar<false>(position,__LINE__);
        assert(currentcc > 0);
        GTTL_SAINUPDATEBUCKETPTR(bucketptr,lastupdatecc,currentcc);
        assert(suftab + suftab_idx < bucketptr);
        *bucketptr++ = (position > 0 &&
                        /* false as used for LONGSEQ */
                        sainseq_getchar<false>(position - 1,__LINE__)
                          < currentcc)
                          ? ~position
                          : position;
#ifdef SAINSHOWSTATE
        assert(bucketptr != nullptr);
        printf("line %d: L-induce: suftab[%zu]=%d\n", __LINE__,
               static_cast<size_t>(bucketptr - 1 - suftab), *(bucketptr - 1));
#endif
      }
    }
  }

  void sain_induceLtypesuffixes2(Sint *suftab)
  {
    if constexpr (T_seqtype == GTTL_SAIN_LONGSEQ)
    {
      sain_LONGSEQ_induceLtypesuffixes2(suftab);
    } else
    {
      sain_PLAIN_MULTISEQ_induceLtypesuffixes2(suftab);
    }
  }

  void sain_special_singleSinduction2(Sint *suftab, Sint position)
  {
    constexpr const bool special_as_pos = T_seqtype == GTTL_SAIN_MULTISEQ;
    const size_t currentcc = sainseq_getchar<special_as_pos>(position,__LINE__);
    if (currentcc < numofchars)
    {
      const size_t putidx = static_cast<size_t>(--bucketfillptr[currentcc]);

      assert(putidx < nonspecialentries);
      suftab[putidx] = (position == 0 ||
                        sainseq_getchar<special_as_pos>(position - 1,__LINE__)
                          > currentcc)
                          ? ~position
                          : position;
#ifdef SAINSHOWSTATE
      printf("end S-induce: suftab[%zu]=%ld\n", putidx,
             static_cast<long>(suftab[putidx]));
#endif
    }
  }

  void sain_LONGSEQ_induceStypesuffixes2(Sint *suftab)
  {
    assert(T_alphasize == 0);
    size_t lastupdatecc = 0;
    SuftabBaseType *fillptr = bucketfillptr;
    Sint *bucketptr = nullptr;

    sain_special_singleSinduction2(suftab, static_cast<Sint>(totallength - 1));
    for (size_t suftab_idx = nonspecialentries; suftab_idx > 0; /* Nothing */)
    {
      suftab_idx--;
      Sint position = suftab[suftab_idx];
      if (position > 0)
      {
        position--;
        /* false as sequence is plainseq */
        const size_t currentcc = sainseq_getchar<false>(position,__LINE__);
        GTTL_SAINUPDATEBUCKETPTR(bucketptr,lastupdatecc,currentcc);
        assert(bucketptr - 1 < suftab + suftab_idx);
        *(--bucketptr) = (position == 0 ||
                          sainseq_getchar<false>(position - 1,__LINE__)
                            > currentcc)
                            ? ~position
                            : position;
#ifdef SAINSHOWSTATE
        assert(bucketptr != nullptr);
        printf("line %d: S-induce: suftab[%zu]=%d\n", __LINE__,
               static_cast<size_t>(bucketptr - suftab), *bucketptr);
#endif
      } else
      {
        suftab[suftab_idx] = ~position;
      }
    }
  }

  void sain_PLAIN_MULTISEQ_induceStypesuffixes2(Sint *suftab)
  {
    assert(T_alphasize == numofchars);
    size_t lastupdatecc = 0;
    SuftabBaseType *fillptr = bucketfillptr;

    sain_special_singleSinduction2(suftab, static_cast<Sint>(totallength - 1));
    if constexpr (T_seqtype == GTTL_SAIN_MULTISEQ)
    {
      induceStypes2fromspecialranges(suftab);
    }
    Sint *bucketptr = nullptr;
    for (size_t suftab_idx = nonspecialentries; suftab_idx > 0; /* Nothing */)
    {
      suftab_idx--;
      Sint position = suftab[suftab_idx];
      if (position > 0)
      {
        position--;
        const size_t currentcc = sainseq_getchar<false>(position,__LINE__);
        /* It seems that currentc is always < T_alphasize */
        GTTL_SAINUPDATEBUCKETPTR(bucketptr,lastupdatecc,currentcc);
        assert(bucketptr - 1 < suftab + suftab_idx);
        constexpr const bool special_as_pos = T_seqtype == GTTL_SAIN_MULTISEQ;
        *(--bucketptr) = (position == 0 ||
                          /* false as sequence is plainseq */
                          sainseq_getchar<special_as_pos>(position - 1,
                                                          __LINE__) > currentcc)
                            ? ~position
                            : position;
#ifdef SAINSHOWSTATE
        printf("line %d: S-induce: suftab[%zu]=%d\n", __LINE__,
               static_cast<size_t>(bucketptr - suftab), *bucketptr);
#endif
      } else
      {
        suftab[suftab_idx] = ~position;
      }
    }
  }

  void sain_induceStypesuffixes2(Sint *suftab)
  {
    if constexpr (T_seqtype == GTTL_SAIN_LONGSEQ)
    {
      sain_LONGSEQ_induceStypesuffixes2(suftab);
    } else
    {
      sain_PLAIN_MULTISEQ_induceStypesuffixes2(suftab);
    }
  }

  void induceStypes2fromspecialranges(Sint *suftab)
  {
    static_assert(T_seqtype == GTTL_SAIN_MULTISEQ);
    assert(T_alphasize == numofchars && seq->sequences_number_get() > 0);

    // get the position of the padding chars and wildcards from second last to
    // first
    constexpr const bool forward = false, invert = true;
    if constexpr (T_alphasize == size_t(4))
    {
      using ThisCharRange = GttlCharRange<char_finder::EncodedNucleotideFinder,
                                          encoded_nucleotide_finder,
                                          forward,
                                          invert>;
      ThisCharRange ranger(multiseq_seq_ptr, totallength);
      for (auto &&range : ranger)
      {
        if (std::get<0>(range) > 0)
        {
          sain_special_singleSinduction2(suftab,
                                         static_cast<Sint>
                                                    (std::get<0>(range)-1));
        }
      }
    } else
    {
      static_assert(T_alphasize == size_t(20));
      using ThisCharRange = GttlCharRange<char_finder::EncodeAminoAcidFinder,
                                          encoded_aa_finder, forward, invert>;
      ThisCharRange ranger(multiseq_seq_ptr, totallength);
      for (auto &&range : ranger)
      {
        if (std::get<0>(range) > 0)
        {
          sain_special_singleSinduction2(suftab,
                                         static_cast<Sint>
                                                    (std::get<0>(range)-1));
        }
      }
    }
  }

  void sain_adjustsuftab(Sint *suftab)
  {
    for (size_t suftab_idx = nonspecialentries; suftab_idx > 0; /* Nothing */)
    {
      suftab_idx--;
      if (suftab[suftab_idx] > 0 &&
          suftab[suftab_idx] < static_cast<Sint>(totallength))
      {
        suftab[suftab_idx] += static_cast<Sint>(totallength);
        Sint *nextgteq = suftab + suftab_idx - 1;

        while (nextgteq >= suftab && *nextgteq < static_cast<Sint>(totallength))
        {
          nextgteq--;
        }
        assert(nextgteq >= suftab);
        if (*nextgteq >= static_cast<Sint>(totallength))
        {
          *nextgteq -= static_cast<Sint>(totallength);
        }
        suftab_idx = static_cast<size_t>(nextgteq - suftab);
      }
    }
  }

  void sain_moveSstar2front(size_t countSstartype, Sint *suftab)
  {
    size_t readidx, writeidx = 0;
    Sint position;

    for (readidx = 0; (position = suftab[readidx]) < 0; readidx++)
    {
      position = ~position;
      suftab[readidx] = position;
      assert(readidx + 1 < nonspecialentries);
    }
    if (readidx < countSstartype)
    {
      for (writeidx = readidx, readidx++; /* Nothing */; readidx++)
      {
        assert(readidx < nonspecialentries);
        if ((position = suftab[readidx]) < 0)
        {
          position = ~position;
          assert(writeidx < readidx);
          suftab[writeidx++] = position;
          suftab[readidx] = 0;

          if (writeidx == countSstartype)
          {
            break;
          }
        } else
        {
          suftab[readidx] = 0;
        }
      }
    } else
    {
#ifndef NDEBUG
      writeidx = readidx;
#endif
    }
    assert(writeidx == countSstartype);
  }

  size_t sain_fast_moveSstar2front(size_t countSstartype, Sint *suftab)
  {
    size_t readidx, namecount = 0, writeidx = 0;
    Sint position;

    for (readidx = 0; (position = suftab[readidx]) < 0; readidx++)
    {
      position = ~position;
      if (position >= static_cast<Sint>(totallength))
      {
        namecount++;
      }
      suftab[readidx] = position;
      assert(readidx + 1 < nonspecialentries);
    }
    if (readidx < countSstartype)
    {
      for (writeidx = readidx, readidx++; /* Nothing */; readidx++)
      {
        assert(readidx < nonspecialentries);
        if ((position = suftab[readidx]) < 0)
        {
          position = ~position;
          if (position >= static_cast<Sint>(totallength))
          {
            namecount++;
          }
          assert(writeidx < readidx);
          suftab[writeidx++] = position;
          suftab[readidx] = 0;
          if (writeidx == countSstartype)
          {
            break;
          }
        } else
        {
          suftab[readidx] = 0;
        }
      }
    } else
    {
#ifndef NDEBUG
      writeidx = readidx;
#endif
    }
    assert(writeidx == countSstartype);
    return namecount;
  }

  void sain_fast_assignSstarnames(size_t countSstartype, SuftabBaseType *suftab,
                                  size_t number_of_names)
  {
    SuftabBaseType *secondhalf = suftab + countSstartype;

    if (number_of_names < countSstartype)
    {
      size_t currentname = number_of_names + 1;
      for (size_t suftab_idx = nonspecialentries; suftab_idx > 0; /* Nothing */)
      {
        suftab_idx--;
        SuftabBaseType position = suftab[suftab_idx];

        if (position >= totallength)
        {
          position -= totallength;
          assert(currentname > 0);
          currentname--;
        }
        if (currentname <= number_of_names)
        {
          secondhalf[position/2] = static_cast<SuftabBaseType>(currentname);
        }
      }
    } else
    {
      for (size_t suftab_idx = 0; suftab_idx < nonspecialentries; suftab_idx++)
      {
        if (suftab[suftab_idx] >= totallength)
        {
          suftab[suftab_idx] -= totallength;
        }
      }
    }
  }

  template<bool fwd>
  void sain_set_undefined(SuftabBaseType *suftab,size_t start, size_t end)
  {
    if constexpr (fwd)
    {
      for (size_t suftab_idx = start; suftab_idx <= end; suftab_idx++)
      {
        suftab[suftab_idx] = 0;
      }
    } else
    {
      assert(start <= end);
      size_t suftab_idx = end;
      while(true)
      {
        suftab[suftab_idx] = 0;
        if (suftab_idx > start)
        {
          suftab_idx--;
        } else
        {
          break;
        }
      }
    }
  }

  void sain_movenames2front(SuftabBaseType *suftab, size_t numberofsuffixes)
  {
    size_t w_idx = numberofsuffixes;
    for (size_t suftab_idx = numberofsuffixes;
         w_idx < 2 * numberofsuffixes; suftab_idx++)
    {
      if (suftab[suftab_idx] > 0)
      {
        suftab[w_idx++] = suftab[suftab_idx] - 1; /* As we have used names with
                                                     offset 1 to distinguish
                                                     them from the undefined
                                                     values signified by 0 */
      }
    }
  }

  void sain_fill_tail_suffixes(SuftabBaseType *suftabtail)
  {
    assert(T_seqtype == GTTL_SAIN_MULTISEQ);
    assert(seq->sequences_number_get() > 0);
    constexpr const bool forward = true,
                         invert = true;
    size_t countspecial = 0;
    if constexpr (T_alphasize == size_t(4))
    {
      using ThisCharRange = GttlCharRange<char_finder::EncodedNucleotideFinder,
                                          encoded_nucleotide_finder,
                                          forward,
                                          invert>;
      ThisCharRange ranger(multiseq_seq_ptr, totallength);
      for (auto &&range : ranger)
      {
        for (size_t idx = 0; idx < std::get<1>(range); idx++)
        {
          suftabtail[countspecial++] = std::get<0>(range) + idx;
        }
      }
    } else
    {
      static_assert(T_alphasize == size_t(20));
      using ThisCharRange = GttlCharRange<char_finder::EncodeAminoAcidFinder,
                                          encoded_aa_finder, forward, invert>;
      ThisCharRange ranger(multiseq_seq_ptr, totallength);
      for (auto &&range : ranger)
      {
        for (size_t idx = 0; idx < std::get<1>(range); idx++)
        {
          suftabtail[countspecial++] = std::get<0>(range) + idx;
        }
      }
    }
  }

  public:
  // Constructor to use for plainseq
  GttlSainseq(GttlMemoryTracker *_memory_tracker,
              const T_seq *_seq, size_t len, const size_t *_charcount_arr)
    : memory_tracker(_memory_tracker)
    , totallength(len)
    , numofchars(T_alphasize)
    , nonspecialentries(len)
    , charcount_arr(_charcount_arr)
    , currentround(0)
    , bucketsizeptr(GTTL_TRACK_MALLOC(SuftabBaseType,
                                      T_alphasize * sizeof *bucketsizeptr))
    , bucketfillptr(GTTL_TRACK_MALLOC(SuftabBaseType,
                                      T_alphasize * sizeof *bucketfillptr))
    , sstarfirstcharcount(GTTL_TRACK_CALLOC(SuftabBaseType,
                                            T_alphasize,
                                            sizeof *sstarfirstcharcount))
    , seq(_seq)
    , bucketsizepoints2suftab(false)
    , bucketfillptrpoints2suftab(false)
    , roundtablepoints2suftab(false)
    , multiseq_seq_ptr(nullptr)
    , this_size(sizeof(GttlSainseq<SuftabBaseType,T_seqtype,T_alphasize>) +
                       T_alphasize * (sizeof *bucketsizeptr +
                                      sizeof *bucketfillptr +
                                      sizeof *sstarfirstcharcount))
  {
    static_assert(T_seqtype == GTTL_SAIN_PLAINSEQ);
    if (use_fast_method())
    {
      roundtable = GTTL_TRACK_MALLOC(SuftabBaseType,
                                     2 * T_alphasize * sizeof *roundtable);
      this_size += 2 * T_alphasize * sizeof *roundtable;
    } else
    {
      roundtable = nullptr;
    }
    for (size_t charidx = 0; charidx < T_alphasize; charidx++)
    {
      bucketsizeptr[charidx] = charcount_arr[charidx];
    }
  }

  // Constructor to use for multiseqs
  GttlSainseq(GttlMemoryTracker *_memory_tracker,
              const T_seq *multiseq, const size_t *_charcount_arr)
    : memory_tracker(_memory_tracker)
    , totallength(multiseq->sequences_total_length_get() +
                  multiseq->sequences_number_get() - 1)
    , numofchars(T_alphasize)
    // the total length of all sequences without the padding_chars and
    // wildcards len including the padding_chars and wildcards
    , nonspecialentries(multiseq->sequences_total_length_get() -
                        _charcount_arr[T_alphasize])
    , charcount_arr(_charcount_arr)
    , currentround(0)
    , bucketsizeptr(GTTL_TRACK_MALLOC(SuftabBaseType,
                                      T_alphasize * sizeof *bucketsizeptr))
    , bucketfillptr(GTTL_TRACK_MALLOC(SuftabBaseType,
                                      T_alphasize * sizeof *bucketfillptr))
    , sstarfirstcharcount(GTTL_TRACK_CALLOC(SuftabBaseType,T_alphasize,
                                            sizeof *sstarfirstcharcount))
    , seq(multiseq)
    , bucketsizepoints2suftab(false)
    , bucketfillptrpoints2suftab(false)
    , roundtablepoints2suftab(false)
    , multiseq_seq_ptr(multiseq->sequence_ptr_get())
    , this_size(sizeof(GttlSainseq<SuftabBaseType,T_seqtype,T_alphasize>) +
                       T_alphasize * (sizeof *bucketsizeptr +
                                      sizeof *bucketfillptr +
                                      sizeof *sstarfirstcharcount))
  {
    assert(_memory_tracker != nullptr);
    static_assert(T_seqtype == GTTL_SAIN_MULTISEQ && T_alphasize > 0);
#ifndef WITH_BRANCHING
    for (size_t charidx = 0; charidx <= T_alphasize; charidx++)
    {
      return_values[charidx] = charidx;
    }
#endif
    if (use_fast_method())
    {
      roundtable = GTTL_TRACK_MALLOC(SuftabBaseType,
                                     2 * T_alphasize * sizeof *roundtable);
      this_size += 2 * T_alphasize * sizeof *roundtable;
    } else
    {
      roundtable = nullptr;
    }
    for (size_t charidx = 0; charidx < T_alphasize; charidx++)
    {
      bucketsizeptr[charidx] = charcount_arr[charidx];
    }
  }

  // Constructor to use for sainseq from array
  GttlSainseq(GttlMemoryTracker *_memory_tracker,
              const T_seq *arr,
              size_t len,
              size_t _numofchars,
              SuftabResources<SuftabBaseType> *suftab_resources_ptr)
    : memory_tracker(_memory_tracker)
    , totallength(len)
    , numofchars(_numofchars)
    , nonspecialentries(len)
    , charcount_arr(nullptr)
    , currentround(0)
    , bucketsizeptr(nullptr)
    , sstarfirstcharcount(nullptr)
    , seq(arr)
    , bucketsizepoints2suftab(false)
    , bucketfillptrpoints2suftab(false)
    , multiseq_seq_ptr(nullptr)
    , this_size(0)
  {
    assert(memory_tracker != nullptr);
    auto vec_subtable
      = suftab_resources_ptr->subtable_assigner(use_fast_method(),numofchars);
    assert(vec_subtable.size() >= 2);
    bucketsizeptr = std::get<0>(vec_subtable[0]);
    bucketsizepoints2suftab = std::get<1>(vec_subtable[0]) == 0;
    bucketfillptr = std::get<0>(vec_subtable[1]);
    bucketfillptrpoints2suftab = std::get<1>(vec_subtable[1]) == 0;
    if (use_fast_method())
    {
      assert(vec_subtable.size() == 3);
      roundtable = std::get<0>(vec_subtable[2]);
      roundtablepoints2suftab = std::get<1>(vec_subtable[2]) == 0;
    } else
    {
      roundtablepoints2suftab = false;
      roundtable = nullptr;
    }

    for (size_t charidx = 0; charidx < numofchars; charidx++)
    {
      bucketsizeptr[charidx] = 0;
    }

    for (size_t pos = 0; pos < totallength; pos++)
    {
      /* false as sequence is array */
      const size_t cc = sainseq_getchar<false>(pos,__LINE__);

      assert(cc < numofchars);
      bucketsizeptr[cc]++;
    }
  }

  // Deconstructor
  ~GttlSainseq(void)
  {
    if (!bucketfillptrpoints2suftab)
    {
      GTTL_UNTRACK_ALLOC(bucketfillptr);
      free(bucketfillptr);
    }
    if (!bucketsizepoints2suftab)
    {
      GTTL_UNTRACK_ALLOC(bucketsizeptr);
      free(bucketsizeptr);
    }
    if (!roundtablepoints2suftab && roundtable != nullptr)
    {
      GTTL_UNTRACK_ALLOC(roundtable);
      free(roundtable);
    }
    if constexpr (T_seqtype != GTTL_SAIN_LONGSEQ)
    {
      GTTL_UNTRACK_ALLOC(sstarfirstcharcount);
      free(sstarfirstcharcount);
    }
  }

  size_t size_in_bytes(void) const noexcept
  {
    return this_size;
  }

  void sain_rec_sortsuffixes(GttlMemoryTracker *memory_tracker,
                             FILE *out_fp,
                             size_t level,
                             RunTimeAtLevel *run_time_at_level,
                             SuftabResources<SuftabBaseType>
                               *suftab_resources_ptr,
                             size_t nonspecialentries,
                             bool intermediatecheck,
                             bool check_suftab,
                             bool buffered)
  {
    SuftabBaseType *suftab = suftab_resources_ptr->suftab_ptr_get();
    RunTimeClass rt_this_level{};
    run_time_at_level->add_level();
    if (out_fp != nullptr)
    {
      fprintf(out_fp,
              "# %zu\t%zu\t%zu\t%.2f\n",
              level, totallength, numofchars,
              static_cast<double>(numofchars) / totallength);
    }
    size_t countSstartype;
    SainbufferKeyValues sain_buffer_key_values(sizeof(SuftabBaseType),
                                               numofchars,totallength);
    if (buffered && sain_buffer_key_values.has_own_memory())
    {
      if (out_fp != nullptr)
      {
        fprintf(out_fp,"# SPACE\tsainbuffer for level %zu (MB):\t%zu\n",
                level,
                static_cast<size_t>(mega_bytes(sain_buffer_key_values
                                                 .size_in_bytes())));
      }
      countSstartype = sain_buffered_insertSstarsuffixes(sain_buffer_key_values,
                                                         suftab);
    } else
    {
      countSstartype = sain_insertSstarsuffixes(suftab);
    }
    if (countSstartype > 0)
    {
      if (roundtable != nullptr)
      {
        sain_incrementfirstSstar(suftab);
      }
      sain_startbuckets();
      sain_induceLtypesuffixes1(reinterpret_cast<Sint *>(suftab));
      if (roundtable != nullptr)
      {
        sain_adjustsuftab(reinterpret_cast<Sint *>(suftab));
      }
      sain_endbuckets();
      sain_induceStypesuffixes1(reinterpret_cast<Sint *>(suftab));
      size_t number_of_names;
      if (roundtable == nullptr)
      {
        sain_moveSstar2front(countSstartype, reinterpret_cast<Sint *>(suftab));
        sain_assignSstarlength(suftab + countSstartype);
        number_of_names = sain_assignSstarnames(countSstartype, suftab);
      } else
      {
        number_of_names = sain_fast_moveSstar2front(countSstartype,
                                                    reinterpret_cast<Sint *>
                                                                    (suftab));
        if (!roundtablepoints2suftab)
        {
          GTTL_UNTRACK_ALLOC(roundtable);
          free(roundtable);
          roundtable = nullptr;
        }
        sain_fast_assignSstarnames(countSstartype, suftab, number_of_names);
      }
      assert(number_of_names <= countSstartype);
      if (number_of_names < countSstartype)
      {
        /* Now the name sequence is in the range from
          countSstartype .. 2 * countSstartype - 1 */
        SuftabBaseType *subseq = suftab + countSstartype;

        /* and we set all what is left of it to undefined */
        sain_set_undefined<true>(suftab, 0, countSstartype - 1);
        if (level == 0)
        {
          /* from the inital call */
          suftab_resources_ptr->firstusable_set(2 * countSstartype);
        }
        sain_movenames2front(suftab, countSstartype);
        GttlSainseq<SuftabBaseType, GTTL_SAIN_LONGSEQ>
          sainseq_rec(memory_tracker,subseq, countSstartype, number_of_names,
                      suftab_resources_ptr);
        memory_tracker->track(&sainseq_rec,__FILE__,__LINE__,
                              sainseq_rec.size_in_bytes());
        if (out_fp != nullptr)
        {
          const size_t mb = static_cast<size_t>(mega_bytes(sainseq_rec
                                                           .size_in_bytes()));
          if (mb > 0)
          {
            fprintf(out_fp,"# SPACE\tlookup tables at level %zu (MB):\t%zu\n",
                    level + 1, mb);
          }
        }
        /* The recursive call at level + 1 */
        sainseq_rec.sain_rec_sortsuffixes(memory_tracker,
                                          out_fp,
                                          level + 1,
                                          run_time_at_level,
                                          suftab_resources_ptr,
                                          countSstartype, /* variable */
                                          intermediatecheck,
                                          check_suftab,
                                          buffered);
        if (out_fp != nullptr)
        {
          run_time_at_level->show(out_fp, level+1);
        }
        sain_expandorder2original(countSstartype, suftab);
        memory_tracker->untrack(&sainseq_rec,__FILE__,__LINE__);
      } else
      {
        if constexpr (T_seqtype == GTTL_SAIN_LONGSEQ)
        {
          assert(sstarfirstcharcount == nullptr);
          sstarfirstcharcount = bucketfillptr;
          for (size_t charidx = 0; charidx < numofchars; charidx++)
          {
            sstarfirstcharcount[charidx] = 0;
            bucketsizeptr[charidx] = 0;
          }
          sain_determineSstarfirstchardist();
        }
      }
    }
    if (intermediatecheck && countSstartype > 0)
    {
      sain_checkorder(suftab, 0, countSstartype - 1);
    }
    if (countSstartype > 0)
    {
      sain_insertsortedSstarsuffixes(suftab, countSstartype - 1);
    }
    sain_startbuckets();
    sain_induceLtypesuffixes2(reinterpret_cast<Sint *>(suftab));
    sain_endbuckets();
    sain_induceStypesuffixes2(reinterpret_cast<Sint *>(suftab));

    if (nonspecialentries > 0 && intermediatecheck)
    {
      sain_checkorder(suftab, 0, nonspecialentries - 1);
    }

    if constexpr (T_seqtype == GTTL_SAIN_MULTISEQ)
    {
      assert(T_alphasize == numofchars);
      sain_fill_tail_suffixes(suftab + nonspecialentries);
      suftab_resources_ptr->place_shortest_suffix();
      if (check_suftab)
      {
        RunTimeClass rt_suftab_check{};
        gttl_suftab_lightweightcheck<SuftabBaseType, char, true>
                                    (seq->sequence_ptr_get(),
                                     totallength,
                                     T_alphasize,
                                     seq->padding_char_get(),
                                     suftab,
                                     numofchars,
                                     charcount_arr);
        rt_suftab_check.show("lightweight check of suftab");
      }
    } else
    {
      if constexpr (T_seqtype == GTTL_SAIN_PLAINSEQ)
      {
        assert(T_alphasize == numofchars);
        suftab_resources_ptr->place_shortest_suffix();
        if (check_suftab)
        {
          RunTimeClass rt_suftab_check{};
          gttl_suftab_lightweightcheck<SuftabBaseType, T_seq, false>
                                      (seq,
                                       totallength,
                                       0, /* wildcard unused */
                                       0, /* paddingchar unused */
                                       suftab,
                                       numofchars,
                                       charcount_arr);
          rt_suftab_check.show("lightweight check of suftab");
        }
      }
    }
    run_time_at_level->set(level,rt_this_level.elapsed());
    if (out_fp != nullptr && level == 0)
    {
      run_time_at_level->show(out_fp, 0);
    }
  }
  size_t nonspecialentries_get(void) const noexcept
  {
    return nonspecialentries;
  }
};

template <typename SuftabBaseType, size_t T_alphasize>
SuftabBaseType *gttl_sain_plain_sortsuffixes(GttlMemoryTracker *memory_tracker,
                                             bool verbose,
                                             const uint8_t *plainseq,
                                             size_t len,
                                             bool intermediatecheck,
                                             bool check_suftab)
{
  assert(memory_tracker != nullptr);
  RunTimeClass rt{};
  size_t *charcount_arr = GTTL_TRACK_CALLOC(size_t,sizeof *charcount_arr,
                                            UCHAR_MAX + 1);
  for (size_t pos = 0; pos < len; pos++)
  {
    charcount_arr[static_cast<int>(plainseq[pos])]++;
  }
  GttlSainseq<SuftabBaseType,GTTL_SAIN_PLAINSEQ,T_alphasize>
             sainseq(memory_tracker,plainseq, len, charcount_arr);
  memory_tracker->track(&sainseq,__FILE__,__LINE__,sainseq.size_in_bytes());
  const size_t level = 0;
  RunTimeAtLevel run_time_at_level;
  SuftabResources<SuftabBaseType> suftab_resources(memory_tracker,len+1);
  if (verbose)
  {
    printf("# SPACE\tsuffixarray (MB):\t%zu\n",
           static_cast<size_t>(mega_bytes(suftab_resources.size_in_bytes())));
    printf("# level\tlen\tasize\tratio\n");
  }
  sainseq.sain_rec_sortsuffixes(memory_tracker,
                                verbose ? stdout : nullptr,
                                level,
                                &run_time_at_level,
                                &suftab_resources,
                                sainseq.nonspecialentries_get(),/*variable*/
                                intermediatecheck,
                                check_suftab,
                                true);
  if (verbose)
  {
    rt.show("construction of suffix array");
  }
  GTTL_UNTRACK_ALLOC(charcount_arr);
  free(charcount_arr);
  memory_tracker->untrack(&sainseq,__FILE__,__LINE__);
  return suftab_resources.suftab_ptr_get();
}

template <typename SuftabBaseType, size_t T_alphasize>
SuftabBaseType *gttl_sain_multiseq_sortsuffixes(GttlMemoryTracker
                                                  *memory_tracker,
                                                bool verbose,
                                                const GttlMultiseq *multiseq,
                                                const size_t *charcount_arr,
                                                bool intermediatecheck,
                                                bool check_suftab,
                                                bool buffered)
{
  assert(memory_tracker != nullptr);
  RunTimeClass rt{};
  GttlSainseq<SuftabBaseType,GTTL_SAIN_MULTISEQ, T_alphasize>
              sainseq(memory_tracker,multiseq,charcount_arr);
  memory_tracker->track(&sainseq,__FILE__,__LINE__,sainseq.size_in_bytes());
  const size_t level = 0;
  RunTimeAtLevel run_time_at_level;
  SuftabResources<SuftabBaseType> suftab_resources
                                    (memory_tracker,
                                     multiseq->sequences_total_length_get() +
                                     multiseq->sequences_number_get() - 1 + 1);
  if (verbose)
  {
    printf("# SPACE\tsuffixarray (MB):\t%zu\n",
           static_cast<size_t>(mega_bytes(suftab_resources.size_in_bytes())));
    const size_t mb = static_cast<size_t>(mega_bytes(sainseq.size_in_bytes()));
    if (mb > 0)
    {
      printf("# SPACE\tlookup tables at level 0 (MB):\t%zu\n",mb);
    }
    printf("# level\tlen\tasize\tratio\n");
  }
  sainseq.sain_rec_sortsuffixes(memory_tracker,
                                verbose ? stdout : nullptr,
                                level,
                                &run_time_at_level,
                                &suftab_resources,
                                sainseq.nonspecialentries_get(),
                                intermediatecheck,
                                check_suftab,
                                buffered);
  if (verbose)
  {
    rt.show("construction of suffix array");
  }
  memory_tracker->untrack(&sainseq,__FILE__,__LINE__);
#ifdef ACCESS_STATISTICS
  size_t total_accesses = 0;
  for (auto &&apl : accesses_per_line)
  {
    printf("line\t%d\t%zu\n",std::get<0>(apl),std::get<1>(apl));
    total_accesses += std::get<1>(apl);
  }
  for (auto &&apl : accesses_special_per_line)
  {
    printf("special line\t%d\t%zu\n",std::get<0>(apl),std::get<1>(apl));
  }
  printf("sequential_count: %zu (%.2f)\n",sequential_count,
                                         static_cast<double>(sequential_count)/
                                         total_accesses);
#endif
  return suftab_resources.suftab_ptr_get();
}
#endif
