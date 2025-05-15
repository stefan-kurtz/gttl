#ifndef QGRAM_DECODER_HPP
#define QGRAM_DECODER_HPP
#include <cstddef>
#include <cassert>
#include "sequences/complement_plain.hpp"

#ifndef NDEBUG
static inline size_t max_integer_code_get(size_t alphabetsize,
                                          size_t qgram_length)
{
  if (alphabetsize == 4)
  {
    if (qgram_length < 32)
    {
      return std::pow(alphabetsize,qgram_length) - 1;
    }
    assert(qgram_length == 32);
    return std::numeric_limits<size_t>::max();
  }
  return gttl_safe_power<size_t>(alphabetsize,qgram_length) - 1;
}
#endif

template<size_t alphabetsize,bool reverse_complement>
class QgramDecoder
{
  private:
  char characters[alphabetsize];
  char *qgram_buffer;
  const size_t qgram_length;
#ifndef NDEBUG
  const size_t max_integer_code;
#endif
  void fill_characters(const char *_characters)
  {
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
  public:
  QgramDecoder(const char *_characters,size_t _qgram_length)
    : qgram_buffer(new char [_qgram_length + 1])
    , qgram_length(_qgram_length)
#ifndef NDEBUG
    , max_integer_code(max_integer_code_get(alphabetsize,_qgram_length))
#endif
  {
    fill_characters(_characters);
  }
  ~QgramDecoder(void)
  {
    delete[] qgram_buffer;
  }
  const char *decode(uint64_t integer_code) const noexcept
  {
    assert(integer_code <= max_integer_code);
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

template <const char *str>
static consteval size_t gttl_strlen(size_t i)
{
  return str[i] == '\0' ? i
                        : gttl_strlen<str>(i+1);
}

template <const char *str>
static consteval size_t gttl_strlen(void)
{
  return gttl_strlen<str>(0);
}

template<const char *characters,bool reverse_complement>
class ConstQgramDecoder
{
  private:
  static constexpr const size_t alphabetsize = gttl_strlen<characters>();
  char *qgram_buffer;
  const size_t qgram_length;
#ifndef NDEBUG
  const size_t max_integer_code;
#endif
  public:
  ConstQgramDecoder(size_t _qgram_length)
    : qgram_buffer(new char [_qgram_length + 1])
    , qgram_length(_qgram_length)
#ifndef NDEBUG
    , max_integer_code(max_integer_code_get(alphabetsize,_qgram_length))
#endif
  {
    static_assert(not reverse_complement or alphabetsize == 4);
  }
  ~ConstQgramDecoder(void)
  {
    delete[] qgram_buffer;
  }
  const char *decode(uint64_t integer_code) const noexcept
  {
    assert(integer_code <= max_integer_code);
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
