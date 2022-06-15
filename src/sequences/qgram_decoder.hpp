#ifndef QGRAM_DECODER_HPP
#define QGRAM_DECODER_HPP
#include <cstddef>
#include <cassert>
#include <cmath>

template<size_t alphabetsize>
class QgramDecoder
{
  private:
  const char *characters;
  size_t qgram_length;
  size_t number_of_all_qgrams;
  char *qgram_buffer;
  public:
  QgramDecoder(const char *_characters,size_t _qgram_length)
    : characters(_characters)
    , qgram_length(_qgram_length)
    , number_of_all_qgrams(std::pow(alphabetsize,qgram_length))
    , qgram_buffer(new char [qgram_length + 1])
    {}
  const char *decode(uint64_t integer_code) const noexcept
  {
    assert(integer_code < number_of_all_qgrams);
    qgram_buffer[qgram_length] = '\0';
    for (char *qgram_ptr = qgram_buffer + qgram_length-1;
         qgram_ptr >= qgram_buffer; qgram_ptr--)
    {
      *qgram_ptr = characters[integer_code % alphabetsize];
      integer_code /= alphabetsize;
    }
    return qgram_buffer;
  };
  ~QgramDecoder(void)
  {
    delete[] qgram_buffer;
  }
};
#endif
