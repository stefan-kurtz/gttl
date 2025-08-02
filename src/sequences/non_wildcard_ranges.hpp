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

#ifndef NON_WILDCARD_RANGES_HPP
#define NON_WILDCARD_RANGES_HPP
#include <cassert>
#include <cstdint>
#include <cstring>
#include <cstddef>
#include <utility>
#include <vector>
#define GTTL_UNALIGNED(X) (size_t(X) & (sizeof (long) - 1))
#define GTTL_UINT64_BYTES  sizeof(uint64_t)

static bool has_null_byte(uint64_t v)
{
  const uint64_t mask = 0x7F7F7F7F7F7F7F7F;
  return  ~((((v & mask) + mask) | v) | mask);
}

inline const void *gttl_memchr(const void *src_void,int search_char,
                               size_t length)
{
  const unsigned char *src = static_cast<const unsigned char *>(src_void);
  const unsigned char delim = static_cast<unsigned char>(search_char);

  while (GTTL_UNALIGNED (src))
  {
    if (!length--)
    {
      return nullptr;
    }
    if (*src == delim)
    {
      return static_cast<const void *>(src);
    }
    src++;
  }

  if (length >= GTTL_UINT64_BYTES)
  {
    /* If we get this far, we know that length is large and src is
       word-aligned. */
    /* The fast code reads the source one word at a time and only
       performs the bytewise search on word-sized segments if they
       contain the search character, which is detected by XORing
       the word-sized segment with a word-sized block of the search
       character and then detecting for the presence of NUL in the
       result.  */
    const uint64_t *asrc = reinterpret_cast<const uint64_t *>(src);
    const uint64_t copied_byte = ~uint64_t(0)/255 * search_char;
    do
    {
      if (has_null_byte(*asrc ^ copied_byte))
      {
        break;
      }
      length -= GTTL_UINT64_BYTES;
      asrc++;
    } while (length >= GTTL_UINT64_BYTES);
    /* If there are fewer than GTTL_UINT64_BYTES characters left,
       then we resort to the bytewise loop.  */
    src = reinterpret_cast<const unsigned char *>(asrc);
  }

  while (length--)
  {
    if (*src == delim)
    {
      return static_cast<const void *>(src);
    }
    src++;
  }
  return nullptr;
}

/* locates the first occurrence in string referenced by src_void, where
   search_char (converted to an unsigned char) does _not_ occur.
   return nullptr if all occurrences in string are search_char.  */

template<char search_char>
static const void *gttl_memcchr(const void *src_void, size_t length)
{
  const unsigned char *src = static_cast<const unsigned char *>(src_void);
  constexpr const unsigned char delim = static_cast<unsigned char>(search_char);

  while (GTTL_UNALIGNED (src))
  {
    if (!length--)
    {
      return nullptr;
    }
    if (*src != delim)
    {
      return static_cast<const void *>(src);
    }
    src++;
  }

  if (length >= GTTL_UINT64_BYTES)
  {
    /* If we get this far, we know that length is large and src is
       word-aligned. */
    /* The fast code reads the source one word at a time and only
       performs the bytewise search on word-sized segments if they
       contain the search character, which is detected by XORing
       the word-sized segment with a word-sized block of the search
       character and then detecting for the presence of NUL in the
       result.  */
    const uint64_t *asrc = reinterpret_cast<const uint64_t *>(src);
    constexpr const uint64_t copied_byte
      = ~uint64_t(0)/uint64_t(255) * static_cast<uint64_t>(search_char);
    do
    {
      if (*asrc != copied_byte)
      {
        break;
      }
      length -= GTTL_UINT64_BYTES;
      asrc++;
    } while (length >= GTTL_UINT64_BYTES);
    /* If there are fewer than GTTL_UINT64_BYTES characters left,
       then we resort to the bytewise loop.  */
    src = reinterpret_cast<const unsigned char *>(asrc);
  }

  while (length--)
  {
    if (*src != delim)
    {
      return static_cast<const void *>(src);
    }
    src++;
  }
  return nullptr;
}

using NonWildCardRangeVector = std::vector<std::pair<size_t,size_t>>;

template<char wildcard>
class NonWildCardRangeIterator
{
  private:
    const char *sequence, *current, *endptr;
    size_t seqlen, minimum_length;
    size_t remaining(const char *curr)
    {
      assert(curr <= endptr);
      return static_cast<size_t>(endptr - curr);
    }
  public:
    NonWildCardRangeIterator(const char *_sequence,size_t _seqlen,
                             size_t _minimum_length = 1)
      : sequence(_sequence),
        current(nullptr),
        endptr(_sequence + _seqlen),
        seqlen(_seqlen),
        minimum_length(_minimum_length)
    {
      assert(_seqlen > 0 && _sequence != nullptr);
    }
    std::vector<std::pair<size_t,size_t>> enumerate(void)
    {
      NonWildCardRangeVector non_wildcard_ranges{};
      if (*sequence == wildcard)
      {
        current = static_cast<const char *>(gttl_memcchr<wildcard>(sequence+1,
                                                                   seqlen - 1));
        if(current == nullptr)
        {
          return non_wildcard_ranges;
        }
      } else
      {
        current = sequence;
      }
      assert(current != nullptr);
      do
      {
        assert(*current != wildcard);
        const char *const next_wildcard = static_cast<const char *>(
                                     memchr(current + 1,
                                            wildcard,
                                            remaining(current + 1)));
        const size_t start = static_cast<size_t>(current - sequence);
        if(next_wildcard == nullptr)
        {
          assert(start <= seqlen - 1);
          if (seqlen - start >= minimum_length)
          {
            non_wildcard_ranges.push_back({start,seqlen-1});
          }
          break;
        }
        const size_t width = static_cast<size_t>(next_wildcard - current);
        assert(width > 0);
        if (width >= minimum_length)
        {
          non_wildcard_ranges.push_back({start,start + width - 1});
        }
        current = static_cast<const char *>(gttl_memcchr<wildcard>
                                             (next_wildcard+1,
                                              remaining(next_wildcard+1)));
      } while(current != nullptr);
      return non_wildcard_ranges;
    }
};
#endif
