#ifndef QGRAM_DECODER_HPP
#define QGRAM_DECODER_HPP
#include <cstddef>
#include <cassert>
#include <cmath>
#include "sequences/complement_plain.hpp"

template<size_t alphabetsize,bool reverse_complement>
class QgramDecoder
{
  private:
  char *characters;
  char *qgram_buffer;
  size_t qgram_length;
#ifndef NDEBUG
  size_t number_of_all_qgrams;
#endif
  public:
  QgramDecoder(const char *_characters,size_t _qgram_length)
    : characters(nullptr)
    , qgram_buffer(new char [_qgram_length + 1])
    , qgram_length(_qgram_length)
    , number_of_all_qgrams(std::pow(alphabetsize,_qgram_length))
  {
    characters = new char [alphabetsize];
    for (size_t idx = 0; idx < alphabetsize; idx++)
    {
      if constexpr (reverse_complement)
      {
        characters[idx] = complement_plain(_characters[idx]);
      } else
      {
        characters[idx] = _characters[idx];
      }
    }
  }
  ~QgramDecoder(void)
  {
    delete[] characters;
    delete[] qgram_buffer;
  }
  const char *decode(uint64_t integer_code) const noexcept
  {
    assert(integer_code < number_of_all_qgrams);
    qgram_buffer[qgram_length] = '\0';
    if constexpr (reverse_complement)
    {
      for (char *qgram_ptr = qgram_buffer;
           qgram_ptr < qgram_buffer + qgram_length; qgram_ptr++)
      {
        *qgram_ptr = characters[integer_code % alphabetsize];
        integer_code /= alphabetsize;
      }
    } else
    {
      for (char *qgram_ptr = qgram_buffer + qgram_length-1;
           qgram_ptr >= qgram_buffer; qgram_ptr--)
      {
        *qgram_ptr = characters[integer_code % alphabetsize];
        integer_code /= alphabetsize;
      }
    }
    return qgram_buffer;
  }
};
#endif
