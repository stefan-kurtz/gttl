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
#ifndef GTTL_LINE_ITERATOR_HPP_
#define GTTL_LINE_ITERATOR_HPP_
#include <cstdbool>
#include <cstdlib>
#include <cassert>
#include <cstring>

#include <string>
#include <algorithm>
#include "utilities/str_format.hpp"
#include "utilities/gttl_file_open.hpp"

template<int buf_size>
class GttlLineIterator
{
  private:
    char buffer[buf_size+1], *bufptr = buffer, *bufend = buffer;
    GttlFpType in_fp;
    char separator;
    bool exhausted,
         endofunit;
    size_t line_number;
    bool fill_buffer(void)
    {
#ifndef QLI_WITHOUT_ZLIB
      const int fill = gzread(in_fp, buffer, static_cast<size_t>(buf_size));
      assert(fill >= 0);
#else
      size_t fill = fread(buffer, 1, static_cast<size_t>(buf_size), in_fp);
#endif

      if (fill > 0)
      {
        bufptr = buffer;
        bufend = buffer + fill;
        assert(fill <= buf_size);
        *bufend = '\0';
        exhausted = static_cast<bool>(fill < buf_size);
        return true;
      }
      exhausted = true;
      return false;
    }
 public:
    GttlLineIterator(GttlFpType _in_fp) :
        in_fp(_in_fp),
        separator(EOF),
        exhausted(false),
        endofunit(false),
        line_number(0)
    {
    }
    void separator_set(char _separator)
    {
      separator = _separator;
    }
    size_t line_number_get(void) const noexcept
    {
      return line_number;
    }
    bool more_lines(void)
    {
      return !(bufptr == bufend && (exhausted || !fill_buffer()));
    }
    bool next(std::string *current_line)
    {
      endofunit = false;
      while (true)
      {
        char *nextnewline, *endptr;
        size_t copy_length;

        if (bufptr == bufend && (exhausted || !fill_buffer()))
        {
          return false;
        }
        nextnewline = reinterpret_cast<char *>
                      (memchr(bufptr, '\n',(size_t) (bufend - bufptr + 1)));
        endptr = nextnewline != nullptr ? nextnewline : (bufend-1);
        copy_length = (size_t) (endptr - bufptr + 1);
        assert(copy_length > 0);
        current_line->append(bufptr, copy_length);
        bufptr = endptr + 1;
        if (nextnewline != nullptr)
        {
          if (bufptr < bufend)
          {
            if (*bufptr == separator)
            {
              endofunit = true;
            }
          } else
          {
            if (!fill_buffer() || *bufptr == separator)
            {
              endofunit = true;
            }
          }
          line_number++;
          return true;
        }
        if (exhausted)
        {
          if (*(bufptr - 1) != '\n')
          {
            StrFormat msg("line %lu: missing newline character",line_number+1);
          }
          endofunit = true;
          line_number++;
          return true;
        }
      }
    }
    bool endofunit_get() const
    {
      return endofunit;
    }
};
#endif
