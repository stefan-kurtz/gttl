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
#ifndef GTTL_LINE_ITERATOR_HPP
#define GTTL_LINE_ITERATOR_HPP
#include <cstdbool>
#include <cstdlib>
#include <cassert>
#include <cstring>

#include <string>
#include <algorithm>
#include <vector>
#include "utilities/str_format.hpp"
#include "utilities/gttl_file_open.hpp"

template<int buf_size>
class GttlLineIterator
{
  private:
    const std::vector<std::string> *inputfiles;
    char buffer[buf_size+1], *bufptr, *bufend;
    GttlFpType in_fp;
    char separator;
    bool file_exhausted,
         endofunit,
         own_in_fp,
         more_files;
    size_t file_index,
           line_number;
    bool fill_buffer(void)
    {
#ifndef QLI_WITHOUT_ZLIB
      const int fill = gzread(in_fp, buffer, static_cast<size_t>(buf_size));
      assert(fill >= 0 && fill <= buf_size);
#else
      size_t fill = fread(buffer, 1, static_cast<size_t>(buf_size), in_fp);
      assert(fill <= static_cast<size_t>(buf_size));
#endif

      if (fill > 0)
      {
        bufptr = buffer;
        bufend = buffer + fill;
        *bufend = '\0';
#ifndef QLI_WITHOUT_ZLIB
        file_exhausted = static_cast<bool>(fill < buf_size);
#else
        file_exhausted = static_cast<bool>(fill <
                                           static_cast<size_t>(buf_size));
#endif
        return true;
      }
      file_exhausted = true;
      return false;
    }
 public:
    GttlLineIterator(GttlFpType _in_fp) :
        inputfiles(nullptr),
        bufptr(buffer),
        bufend(buffer),
        in_fp(_in_fp),
        separator(EOF),
        file_exhausted(false),
        endofunit(false),
        own_in_fp(false),
        more_files(false),
        file_index(0),
        line_number(0)
    {
    }
    GttlLineIterator(const char *inputfile) :
        inputfiles(nullptr),
        bufptr(buffer),
        bufend(buffer),
        separator(EOF),
        file_exhausted(false),
        endofunit(false),
        own_in_fp(true),
        more_files(false),
        file_index(0),
        line_number(0)
    {
      in_fp = gttl_fp_type_open(inputfile,"rb");
      if (in_fp == nullptr)
      {
        throw std::string(": cannot open file");
      }
    }
    GttlLineIterator(const std::vector<std::string> *_inputfiles) :
        inputfiles(_inputfiles),
        bufptr(buffer),
        bufend(buffer),
        separator(EOF),
        file_exhausted(false),
        endofunit(false),
        own_in_fp(true),
        file_index(0),
        line_number(0)
    {
      assert(_inputfiles->size() > 0);
      more_files = inputfiles->size() > 1;
      in_fp = gttl_fp_type_open(inputfiles->at(0).c_str(),"rb");
      if (in_fp == nullptr)
      {
        throw std::string(": cannot open file");
      }
    }
    ~GttlLineIterator(void)
    {
      if (own_in_fp)
      {
        gttl_fp_type_close(in_fp);
      }
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
      if (bufptr != bufend || (!file_exhausted && fill_buffer()))
      {
        return true;
      }
      if (more_files)
      {
        bufptr = bufend = buffer;
        file_exhausted = false;
        endofunit = false;
        gttl_fp_type_close(in_fp);
        line_number = 0;
        file_index++;
        assert(file_index < inputfiles->size());
        more_files = file_index + 1 < inputfiles->size();
        in_fp = gttl_fp_type_open(inputfiles->at(file_index).c_str(),"rb");
        if (in_fp == nullptr)
        {
          throw std::string(": cannot open file");
        }
        return true;
      }
      return false;
    }
    bool next(std::string *current_line)
    {
      endofunit = false;
      while (more_lines())
      {
        char *nextnewline = reinterpret_cast<char *>
                            (memchr(bufptr,'\n',
                                    static_cast<size_t>(bufend - bufptr + 1)));
        if (nextnewline != nullptr)
        {
          const size_t copy_length = static_cast<size_t>(nextnewline + 1
                                                         - bufptr);
          assert(copy_length > 0);
          current_line->append(bufptr, copy_length);
          bufptr = nextnewline + 1;
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
        const size_t copy_length = static_cast<size_t>(bufend - bufptr);
        if (copy_length > 0)
        {
          current_line->append(bufptr, copy_length);
        }
        bufptr = bufend;
        if (file_exhausted)
        {
          if (*(bufptr - 1) != '\n')
          {
            StrFormat msg(", line %lu: missing newline character",
                          line_number+1);
            throw msg.str();
          }
          endofunit = true;
          line_number++;
          return true;
        }
      }
      return false;
    }
    bool endofunit_get(void) const noexcept
    {
      return endofunit;
    }
    size_t file_index_get(void) const noexcept
    {
      return file_index;
    }
};
#endif
