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
#ifndef GTTL_FASTQ_ITERATOR_HPP
#define GTTL_FASTQ_ITERATOR_HPP
#include <array>
#include "utilities/gttl_line_iterator.hpp"

template <int buf_size>
class GttlFastQIterator
{
  private:
    std::array<std::string,4> seqbufs{};
    int state = 0;
    GttlLineIterator<buf_size> gttl_li;

  public:
    GttlFastQIterator(QliFpType _in_fp) :
      gttl_li(GttlLineIterator<buf_size>(_in_fp))
    {
    }
    bool next(void)
    {
      for (size_t idx = 0; idx < seqbufs.size(); idx++)
      {
        seqbufs[idx].clear();
      }
      while (gttl_li.next(&seqbufs[state]))
      {
        seqbufs[state].pop_back();
        if (state == 3)
        {
          state = 0;
          return true;
        }
        state++;
      }
      return false;
    }

    const std::string *header_get() const
    {
      return &seqbufs[0];
    }

    const std::string *sequence_get() const
    {
      return &seqbufs[1];
    }

    const std::string *quality_get() const
    {
      return &seqbufs[3];
    }
};
#endif
