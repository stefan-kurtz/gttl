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
#include <vector>
#define GTTL_UNALIGNED(X) ((size_t)(X) & (sizeof (long) - 1))
#define GTTL_ULONG_BYTES  sizeof(unsigned long)

static bool has_null_byte(unsigned long v)
{
  const unsigned long mask = 0x7F7F7F7F7F7F7F7F;
  return  ~((((v & mask) + mask) | v) | mask);
}

inline void *gttl_memchr(const void *src_void,int search_char, size_t length)
{
  const unsigned char *src = (const unsigned char *) src_void;
  const unsigned char delim = (char) search_char;

  while (GTTL_UNALIGNED (src))
  {
    if (!length--)
    {
      return nullptr;
    }
    if (*src == delim)
    {
      return (void *) src;
    }
    src++;
  }

  if (length >= GTTL_ULONG_BYTES)
  {
    /* If we get this far, we know that length is large and src is
       word-aligned. */
    /* The fast code reads the source one word at a time and only
       performs the bytewise search on word-sized segments if they
       contain the search character, which is detected by XORing
       the word-sized segment with a word-sized block of the search
       character and then detecting for the presence of NUL in the
       result.  */
    unsigned long *asrc = (unsigned long *) src;
    const unsigned long copied_byte = ~0UL/255 * search_char;
    do
    {
      if (has_null_byte(*asrc ^ copied_byte))
      {
        break;
      }
      length -= GTTL_ULONG_BYTES;
      asrc++;
    } while (length >= GTTL_ULONG_BYTES);
    /* If there are fewer than GTTL_ULONG_BYTES characters left,
       then we resort to the bytewise loop.  */
    src = (unsigned char *) asrc;
  }

  while (length--)
  {
    if (*src == delim)
    {
      return (void *) src;
    }
    src++;
  }
  return nullptr;
}

/* locates the first occurrence in string referenced by src_void, where
   search_char (converted to an unsigned char) does not occur.
   return nullptr if all occurrences in string are search_char. */

static void *gttl_memcchr(const void *src_void,int search_char, size_t length)
{
  const unsigned char *src = (const unsigned char *) src_void;
  const unsigned char delim = (char) search_char;

  while (GTTL_UNALIGNED (src))
  {
    if (!length--)
    {
      return nullptr;
    }
    if (*src != delim)
    {
      return (void *) src;
    }
    src++;
  }

  if (length >= GTTL_ULONG_BYTES)
  {
    /* If we get this far, we know that length is large and src is
       word-aligned. */
    /* The fast code reads the source one word at a time and only
       performs the bytewise search on word-sized segments if they
       contain the search character, which is detected by XORing
       the word-sized segment with a word-sized block of the search
       character and then detecting for the presence of NUL in the
       result.  */
    unsigned long *asrc = (unsigned long *) src;
    const unsigned long copied_byte = ~0UL/255 * search_char;
    do
    {
      if (*asrc != copied_byte)
      {
        break;
      }
      length -= GTTL_ULONG_BYTES;
      asrc++;
    } while (length >= GTTL_ULONG_BYTES);
    /* If there are fewer than GTTL_ULONG_BYTES characters left,
       then we resort to the bytewise loop.  */
    src = (unsigned char *) asrc;
  }

  while (length--)
  {
    if (*src != delim)
    {
      return (void *) src;
    }
    src++;
  }
  return nullptr;
}

typedef std::vector<std::pair<size_t,size_t>> NonWildCardRangeVector;

template<char wildcard>
class NonWildCardRangeIterator
{
  private:
    const char *sequence, *current, *endptr;
    size_t seqlen;
    size_t remaining(const char *curr)
    {
      assert(curr <= endptr);
      return (size_t) (endptr - curr);
    }
  public:
    NonWildCardRangeIterator(const char *_sequence,size_t _seqlen) :
      sequence(_sequence),
      current(nullptr),
      endptr(_sequence + _seqlen),
      seqlen(_seqlen)
    {
      assert(_seqlen > 0 && _sequence != NULL);
    }
    std::vector<std::pair<size_t,size_t>> enumerate(void)
    {
      NonWildCardRangeVector non_wildcard_ranges{};
      if (*sequence == wildcard)
      {
        current = (const char *) gttl_memcchr(sequence+1,wildcard,seqlen - 1);
        if (current == NULL)
        {
          return non_wildcard_ranges;
        }
      } else
      {
        current = sequence;
      }
      assert(current != NULL);
      do
      {
        assert(*current != wildcard);
        const char *next_wildcard
          = (const char *) memchr(current+1,wildcard,remaining(current+1));
        size_t start = (size_t) (current - sequence);
        if (next_wildcard == NULL)
        {
          assert(start <= seqlen - 1);
          non_wildcard_ranges.push_back({start,seqlen-1});
          break;
        }
        size_t width = (size_t) (next_wildcard - current);
        assert(width > 0);
        non_wildcard_ranges.push_back({start,start + width - 1});
        current = (const char *) gttl_memcchr(next_wildcard+1,wildcard,
                                         remaining(next_wildcard+1));
      } while (current != NULL);
      return non_wildcard_ranges;
    }
};
#endif
