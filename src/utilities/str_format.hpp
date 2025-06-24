#ifndef STR_FORMAT_HPP
#define STR_FORMAT_HPP
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
#include <cstdarg>
#include <cstdio>
#include <string>

#define STR_FORMAT_VARIABLE(OPERATOR)\
        va_list vl;\
        va_start(vl, fmt);\
        int size = vsnprintf(0, 0, fmt, vl) + sizeof('\0');\
        va_end(vl);\
        char *buffer = new char [size];\
        va_start(vl, fmt);\
        size = vsnprintf(buffer, size, fmt, vl);\
        va_end(vl);\
        this_string OPERATOR std::string(buffer, size);\
        delete[] buffer

class StrFormat
{
  private:
    std::string this_string{};
  public:
    StrFormat(const char* fmt, ...) __attribute__ ((format (printf, 2, 3)))
    {
      STR_FORMAT_VARIABLE(=);
    }
    std::string str(void) const
    {
      return this_string;
    }
    void append(std::string &more)
    {
      this_string += more;
    }
    void append(const char* fmt, ...) __attribute__ ((format (printf, 2, 3)))
    {
      STR_FORMAT_VARIABLE(+=);
    }
};

#endif
