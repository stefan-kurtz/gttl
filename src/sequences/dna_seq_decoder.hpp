#ifndef DNA_SEQ_DECODER_HPP
#define DNA_SEQ_DECODER_HPP

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <string>
class DNAQgramDecoder
{
  class Iterator
  {
    const uint64_t mask,
                   *sub_unit_ptr;
    size_t current_qgram_idx, idx_of_unit, shift_last;
    uint64_t integer;
    public:
    Iterator(size_t qgram_length,
             size_t _current_qgram_idx,
             const uint64_t *_sub_unit_ptr)
      : mask(((uint64_t(1) << (2 * (qgram_length - 1))) - 1) << 2)
      , sub_unit_ptr(_sub_unit_ptr)
      , current_qgram_idx(_current_qgram_idx)
      , idx_of_unit(_current_qgram_idx + qgram_length < size_t(32) ? 0
                                                                   : size_t(1))
      , shift_last(qgram_length == size_t(32) ? size_t(62)
                                              : (size_t(62) -
                                                 (2 * qgram_length)))
      , integer(sub_unit_ptr[0] >> (size_t(64) - 2 * qgram_length))
    { }
    uint64_t operator*() const
    {
      return integer;
    }
    Iterator& operator++() /* prefix increment*/
    {
      integer = (integer << 2) & mask;
      assert((integer & uint64_t(3)) == 0);
      const uint64_t new_char = (sub_unit_ptr[idx_of_unit] >> shift_last)
                                & uint64_t(3);
      integer |= new_char;
      if (shift_last > 0)
      {
        assert(shift_last >= size_t(2));
        shift_last -= size_t(2);
      } else
      {
        shift_last = size_t(62);
        idx_of_unit++;
      }
      current_qgram_idx++;
      return *this;
    }
    bool operator != (const Iterator& other) const
    {
      return current_qgram_idx != other.current_qgram_idx;
    }
  };
  const size_t qgram_length;
  const size_t number_of_qgrams;
  const uint64_t *sub_unit_ptr;
  public:
  DNAQgramDecoder(size_t _qgram_length,
                  size_t _number_of_qgrams,
                  const uint64_t *_sub_unit_ptr)
    : qgram_length(_qgram_length)
    , number_of_qgrams(_number_of_qgrams)
    , sub_unit_ptr(_sub_unit_ptr)
  {}
  [[nodiscard]] Iterator begin(void) const
  {
    return {qgram_length, 0, sub_unit_ptr};
  }
  [[nodiscard]] Iterator end(void) const
  {
    return {qgram_length, number_of_qgrams, sub_unit_ptr};
  }
};

static inline std::basic_string<uint8_t> dna_sequence_decode(
                                                const uint64_t *sub_unit_ptr,
                                                size_t sequence_length)
{
  int shift = 62;
  size_t current_unit = 0;
  std::basic_string<uint8_t> s;
  s.reserve(sequence_length);
  for (size_t idx = 0; idx < sequence_length; idx++)
  {
    s.push_back(static_cast<uint8_t>(sub_unit_ptr[current_unit] >> shift)
                & uint8_t(3));
    if (shift >= 2)
    {
      shift -= 2;
    } else
    {
      shift = 62;
      current_unit++;
    }
  }
  return s;
}
#endif
