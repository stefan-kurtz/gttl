#ifndef GTTL_LINE_GENERATOR_HPP
#define GTTL_LINE_GENERATOR_HPP

#include <algorithm>
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
    , input_string(nullptr)
  {
    out = default_buffer;
  }

  GttlLineGenerator(const char* file_name, bool _is_end = false)
    : file(gttl_fp_type_open(file_name, "rb"))
    , is_end(_is_end)
    , line_number(0)
    , input_string(nullptr)
  {
    out = default_buffer;
  }

  GttlLineGenerator(const char* file_name,
                    char (&_out)[buf_size],
                    bool _is_end = false)
    : file(gttl_fp_type_open(file_name, "rb"))
    , out(_out)
    , is_end(_is_end)
    , line_number(0)
    , input_string(nullptr)
    {}


  GttlLineGenerator(GttlFpType fp, char (&_out)[buf_size], bool _is_end = false)
    : file(fp)
    , out(_out)
    , is_end(_is_end)
    , line_number(0)
    , input_string(nullptr)
    {}

  GttlLineGenerator(const char* _input_string, size_t _string_length)
    : line_number(0)
    , input_string(_input_string)
    , string_length(_string_length)
    , current_ptr(_input_string)
    {}

  ~GttlLineGenerator()
  {
    gttl_fp_type_close(file);
  }

  bool advance(size_t *length = nullptr)
  {
    if(is_end) return false;
    ++line_number;

    if(input_string != nullptr)
    {
      if(current_ptr >= input_string + string_length)
      {
        is_end = true;
        return false;
      }


      const char* next_newline = std::find(current_ptr,
                                           input_string + string_length,
                                           '\n');
      const char* line_end = next_newline;

      if(out != nullptr)
      {
        const size_t line_len = line_end - current_ptr;
        const size_t copy_len = std::min(line_len, buf_size - 1);
        std::memcpy(out, current_ptr, copy_len);
        out[copy_len] = '\0';
      }

      current_ptr = (next_newline < input_string + string_length)
        ? next_newline + 1
        : input_string + string_length;
      return true;
    }

    if(out == nullptr)
    {
      char to_discard[buf_size];
      if(!gttl_fp_type_gets(file, to_discard, buf_size))
        is_end = true;
      return !is_end;
    }
    if(!gttl_fp_type_gets(file, out, buf_size))
      is_end = true;
    size_t len = std::strlen(out);
    while(len > 0 && (out[len-1] == '\r' || out[len-1] == '\n'))
    {
      out[--len] = '\0';
    }

    if(length != nullptr)
    {
      *length = len;
    }
    return !is_end;
  }

  void reset(void)
  {
    is_end = false;
    line_number = 0;
    if(input_string != nullptr)
    {
      current_ptr = input_string;
    }
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

  const char* operator*() const
  {
    return out;
  }

  Iterator begin()
  {
    return Iterator(this, false);
  }

  Iterator end()
  {
    return Iterator(this, true);
  }

  char getc(void)
  {
    if(input_string != nullptr)
    {
      char ret = *input_string;
      input_string += sizeof(char);
      return ret;
    }
    return gttl_fp_type_getc(file);
  }

  private:
  GttlFpType file;
  char* out;
  bool is_end;
  char default_buffer[buf_size];
  size_t line_number;

  const char* input_string;
  size_t string_length;
  const char* current_ptr;
};

#endif  // GTTL_LINE_GENERATOR_HPP
