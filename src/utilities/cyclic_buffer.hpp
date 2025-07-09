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
#ifndef CYCLIC_BUFFER_HPP
#define CYCLIC_BUFFER_HPP
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cassert>
template<typename Basetype,size_t max_size>
class CyclicBuffer
{
  private:
  size_t max_num_elems, num_elems;
  uint8_t first, next[max_size];
  Basetype space[max_size + max_size - 1]; /* extra space to copy values
                                              from index 0..first-1 to
                                              the end of the buffer,
                                              as required for  */
  public:
  CyclicBuffer(void) :
    max_num_elems(0),
    num_elems(0),
    first(0)
  {
  }
  void initialize(size_t _max_num_elems) noexcept
  {
    max_num_elems = _max_num_elems;
    assert(max_num_elems >= 1 && max_num_elems <= max_size);
    for (size_t idx = 0; idx < max_num_elems; idx++)
    {
      next[idx] = idx < max_num_elems - 1 ? idx + 1 : 0;
    }
  }
  void prepend(Basetype new_elem) noexcept
  {
    assert(num_elems < max_num_elems);
    num_elems++;
    space[max_num_elems - num_elems] = new_elem;
  }
  void append(Basetype new_elem) noexcept
  {
    assert(num_elems < max_num_elems);
    space[num_elems++] = new_elem;
  }
  Basetype shift(Basetype new_elem) noexcept
  {
    assert(max_num_elems == num_elems);
    Basetype ret = space[first];
    space[first] = new_elem;
    first = next[first];
    return ret;
  }
  const Basetype *consecutive_memory_buffer_content(void) noexcept
  {
    if (first > 0)
    {
      /* copy initial elements of space from index 0 to first-1
         after the end of the space used for the cyclic buffer */
      assert(first < max_size);
      memcpy((void *) &space[max_num_elems],(void *) &space[0],
             static_cast<size_t>(first) * sizeof *space);
    }
    return space + static_cast<size_t>(first);
  }
  [[nodiscard]] const Basetype *pointer_to_array(void) const noexcept
  {
    assert(num_elems == max_num_elems);
    return &space[0];
  }
  [[nodiscard]] size_t size(void) const noexcept { return num_elems; }
  void reset(void) noexcept
  {
    num_elems = 0;
    first = 0;
  }
};
#endif
