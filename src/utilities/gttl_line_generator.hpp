#ifndef GTTL_LINE_GENERATOR_HPP
#define GTTL_LINE_GENERATOR_HPP

#include <cstring>
#include "utilities/gttl_file_open.hpp"


template <const size_t buf_size = (1 << 14)>
class GttlLineGenerator
{
  public:
  GttlLineGenerator(GttlFpType fp, bool _is_end = false)
    : file(fp)
    , is_end(_is_end)
    , line_number(0)
  {
    out = default_buffer;
  }

  GttlLineGenerator(const char* file_name, bool _is_end = false)
    : file(gttl_fp_type_open(file_name, "rb"))
    , is_end(_is_end)
    , line_number(0)
  {
    out = default_buffer;
  }

  GttlLineGenerator(const char* file_name,
                    char (&_out)[buf_size],
                    bool _is_end = false)
    : file(gttl_fp_type_open(file_name, "rb"))
    , out(_out)
    , is_end(_is_end)
    , line_number(0) {}


  GttlLineGenerator(GttlFpType fp, char (&_out)[buf_size], bool _is_end = false)
    : file(fp)
    , out(_out)
    , is_end(_is_end)
    , line_number(0) {}

  bool advance(void)
  {
    if(is_end) return false;
    ++line_number;
    if(out == nullptr)
    {
      char to_discard[buf_size];
      gttl_fp_type_gets(file, to_discard, buf_size);
      return true;
    }
    if(!gttl_fp_type_gets(file, out, buf_size))
      is_end = true;
    return !is_end;
  }

  void reset(void)
  {
#ifndef GTTL_WITHOUT_ZLIB
    gzrewind(file);
#else
    rewind(file);
#endif
  }

  size_t line_number_get(void) const noexcept
  {
    return line_number;
  }

  void set_out_buffer(char* _out)
  {
    out = _out;
  }

  class Iterator
  {
    public:
    Iterator(GttlLineGenerator* generator, bool end = false)
      : gen(generator), is_end(end)
    {
      if(!end) ++(*this);
    }

    const char* operator*() const
    {
      return gen->out;
    }

    Iterator& operator++()
    {
      if(!gen->advance())
        is_end = true;
      return *this;
    }

    bool operator==(const Iterator& other) const
    {
      return (is_end == other.is_end) && (gen == other.gen);
    }

    bool operator!=(const Iterator& other) const
    {
      return !(*this == other);
    }

    private:
    GttlLineGenerator* gen;
    bool is_end;
  };

  Iterator begin()
  {
    return Iterator(this, false);
  }

  Iterator end()
  {
    return Iterator(this, true);
  }

  private:
  GttlFpType file;
  char* out;
  bool is_end;
  char default_buffer[buf_size];
  size_t line_number;
};

#endif  // GTTL_LINE_GENERATOR_HPP
