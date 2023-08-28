#ifndef GTTL_SUBSTRING_HPP
#define GTTL_SUBSTRING_HPP

#include <cstddef>
#include <cstdbool>
#include <cassert>
#include "sequences/complement_uint8.hpp"

template<typename CharType>
class GttlSubstring
{
  bool forward_strand;
  const CharType *original_seq;
  size_t start,
         substringlen,
         original_seqlen;
  public:
  GttlSubstring(void)
    : forward_strand(true)
    , original_seq(nullptr)
    , start(0)
    , substringlen(0)
    , original_seqlen(0)
  {}
  GttlSubstring(const CharType *_original_seq,
                size_t _start,
                size_t _substringlen)
    : forward_strand(true)
    , original_seq(_original_seq)
    , start(_start)
    , substringlen(_substringlen)
    , original_seqlen(0)
  {}
  GttlSubstring(bool _forward_strand,
                const CharType *_original_seq,
                size_t _start,
                size_t _substringlen,
                size_t _original_seqlen)
    : forward_strand(_forward_strand)
    , original_seq(_original_seq)
    , start(_start)
    , substringlen(_substringlen)
    , original_seqlen(_original_seqlen)
  {}
  void set(const CharType *_original_seq,
           size_t _start,
           size_t _substringlen)
  {
    forward_strand = true;
    original_seq = _original_seq;
    start = _start;
    substringlen = _substringlen;
    original_seqlen = 0;
  }
  void set(bool _forward_strand,
           const CharType *_original_seq,
           size_t _start,
           size_t _substringlen,
           size_t _original_seqlen)
  {
    forward_strand = _forward_strand;
    original_seq = _original_seq;
    start = _start;
    substringlen = _substringlen;
    original_seqlen = _original_seqlen;
  }
  CharType operator [](size_t idx) const noexcept
  {
    assert(idx < substringlen);
    if (forward_strand)
    {
      return original_seq[start + idx];
    }
    assert(start < original_seqlen);
    const size_t transformed_end = original_seqlen - 1 - start;
    assert(idx <= transformed_end);
    return complement_uint8(original_seq[transformed_end - idx]);
  }
  std::string to_string(void) const noexcept
  {
    if (forward_strand)
    {
      return "+," + std::to_string(start) + "," + std::to_string(substringlen);
    }
    return "-," + std::to_string(start) + ","
                + std::to_string(substringlen) + ","
                + std::to_string(original_seqlen);
  }
  size_t size(void) const noexcept
  {
    return substringlen;
  }
};
#endif
