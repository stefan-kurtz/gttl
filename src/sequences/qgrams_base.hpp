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
#ifndef QGRAMS_HPP
#define QGRAMS_HPP
#include <cstddef>
#include <cstdint>
#include <cmath>
void qgram_show(const uint8_t *qgram,size_t qgram_length);

void qgram_decode(size_t alpha_size,uint8_t *qgram,size_t qgram_length,
                  uint64_t integer_code);

uint64_t qgram_encode(size_t alpha_size,const uint8_t *qgram,
                      size_t qgram_length);

template<typename T,void (*process)(uint64_t,T *)>
void qgram_enumerate_callback(size_t alpha_size,
                              size_t qgram_length,
                              const uint8_t *sequence,
                              size_t seqlen,
                              T *data)
{
  if (qgram_length <= seqlen)
  {
    const uint64_t multiplier = (uint64_t) std::pow(alpha_size,qgram_length-1);
    uint64_t integer_code = qgram_encode(alpha_size,sequence,qgram_length);
    process(integer_code,data);
    for (size_t idx = 0; idx < seqlen - qgram_length; idx++)
    {
      integer_code -= multiplier * sequence[idx];
      integer_code *= alpha_size;
      integer_code += sequence[idx+qgram_length];
      process(integer_code,data);
    }
  }
}
#endif
