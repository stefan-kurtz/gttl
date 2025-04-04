#ifndef GTTL_LINE_GENERATOR_HPP
#define GTTL_LINE_GENERATOR_HPP

#include "utilities/gttl_file_open.hpp"
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <string>
#include <type_traits>

template <const size_t buf_size = (1U << 14U), const bool use_heap = false>
class GttlLineGenerator
{
  public:
  explicit GttlLineGenerator (GttlFpType fp, bool _is_end = false)
    : file(fp)
    , is_end(_is_end)
    , line_number(0)
  {
    set_default_out_buffer();
  }

  explicit GttlLineGenerator(const char* file_name, bool _is_end = false)
    : file(gttl_fp_type_open(file_name, "rb"))
    , is_end(_is_end)
    , line_number(0)
  {
    set_default_out_buffer();
  }

  explicit GttlLineGenerator(const char* file_name,
                    std::string* _out,
                    bool _is_end = false) requires use_heap
    : file(gttl_fp_type_open(file_name, "rb"))
    , out(_out)
    , is_end(_is_end)
    , line_number(0)
    {}

  explicit GttlLineGenerator(const char* file_name,
                    char (&_out)[buf_size],
                    bool _is_end = false) requires (not use_heap)
    : file(gttl_fp_type_open(file_name, "rb"))
    , out(_out)
    , is_end(_is_end)
    , line_number(0)
    {}

  explicit GttlLineGenerator(GttlFpType fp,
                             std::string* _out,
                             bool _is_end = false) requires use_heap
    : file(fp)
    , out(_out)
    , is_end(_is_end)
    , line_number(0)
    {}

  explicit GttlLineGenerator(GttlFpType fp,
                             char (&_out)[buf_size],
                             bool _is_end = false) requires (not use_heap)
    : file(fp)
    , out(_out)
    , is_end(_is_end)
    , line_number(0)
    {}

  explicit GttlLineGenerator(const char* _input_string, size_t _string_length)
    : file(nullptr)
    , is_end(_string_length == 0)
    , line_number(0)
    , input_string(_input_string)
    , string_length(_string_length)
    , current_ptr(_input_string)
    {
      set_default_out_buffer();
    }

  // Delete copy/move constructur & assignment operator
  GttlLineGenerator(const GttlLineGenerator &other) = delete;
  GttlLineGenerator& operator=(const GttlLineGenerator &other) = delete;
  GttlLineGenerator(GttlLineGenerator&& other) = delete;
  GttlLineGenerator& operator=(GttlLineGenerator&& other) = delete;

  ~GttlLineGenerator()
  {
    gttl_fp_type_close(file);
  }

private:
  bool read_from_mapped_string(size_t *length)
  {
    if(current_ptr >= input_string + string_length)
    {
      is_end = true;
      return false;
    }

    const char* next_newline = std::find(current_ptr,
                                         input_string + string_length,
                                         '\n');
    const size_t line_len = next_newline - current_ptr;
    const size_t copy_len = std::min(line_len, buf_size - 1);

    if(out != nullptr)
    {
      if constexpr (use_heap)
      {
        out->resize(copy_len);
        std::memcpy(out->data(), current_ptr, copy_len);
      }else
      {
        //TODO: Don't use 0-terminated string here
        std::memcpy(out, current_ptr, copy_len);
        out[copy_len] = '\0';
      }
    }

    if(length != nullptr)
    {
      (*length) = copy_len;
    }

    current_ptr = (next_newline < input_string + string_length)
      ? next_newline + 1
      : input_string + string_length;
    return true;
  }

  bool read_from_file(size_t* length)
  {
    // ch is an int, which is larger than char, to accomodate
    // error-values and EOF (=-1)
    int ch = EOF;
    size_t len = 0;
    while((ch = gttl_fp_type_getc(file)) != EOF)
    {
      if(ch == '\r') continue;
      if(ch == '\n') break;
      if constexpr(use_heap)
      {
        out->push_back(static_cast<char>(ch));
      }else
      {
        out[len] = static_cast<char>(ch);
        len++;
      }
    }
    if constexpr(not use_heap)
    {
      out[len] = '\0';
    }
    if(ch == EOF) is_end = true;

    if(length != nullptr)
    {
      if constexpr(use_heap)
      {
        *length = out->size();
      }else
      {
        *length = len;
      }
    }
    return not is_end;
  }

  bool discard_line(size_t *length = nullptr)
  {
    char ch = EOF;
    if(length != nullptr)
    {
      *length = 0;
    }
    while((ch = gttl_fp_type_getc(file)) != EOF)
    {
      if(ch == '\n') break;
      if(length != nullptr)
      {
        (*length)++;
      }
    }
    return (ch != EOF);
  }

public:
  bool advance(size_t *length = nullptr)
  {
    if(length != nullptr)
    {
      *length = 0;
    }
    if(is_end) return false;
    ++line_number;

    if constexpr(use_heap)
    {
      if(out != nullptr) out->clear();
    }

    if(input_string != nullptr)
    {
      return read_from_mapped_string(length);
    }

    if(out == nullptr)
    {
      return discard_line(length);
    }

    return read_from_file(length);
  }

  void reset()
  {
    is_end = false;
    line_number = 0;
    if(input_string != nullptr)
    {
      current_ptr = input_string;
    }
    gttl_fp_type_rewind(file);
  }

  [[nodiscard]] size_t line_number_get() const noexcept
  {
    return line_number;
  }

  void set_out_buffer(std::string* _out) requires use_heap
  {
    out = _out;
  }

  void set_out_buffer(char* _out) requires (not use_heap)
  {
    out = _out;
  }

  void set_default_out_buffer()
  {
    if constexpr(use_heap)
    {
      out = &default_buffer;
    }else
    {
      out = default_buffer;
    }
  }

  class Iterator
  {
    public:
    explicit Iterator(GttlLineGenerator* generator, bool end = false)
      : gen(generator), is_end(end)
    {
      if(not end) ++(*this);
    }

    const std::string& operator*() const requires use_heap
    {
      return *(gen->out);
    }

    const char* operator*() const requires (not use_heap)
    {
      return gen->out;
    }

    Iterator& operator++()
    {
      is_end = not gen->advance();
      return *this;
    }

    bool operator==(const Iterator& other) const
    {
      return (is_end == other.is_end) and (gen == other.gen);
    }

    private:
    GttlLineGenerator* gen;
    bool is_end;
  };

  std::string operator*() const requires use_heap
  {
    return *out;
  }

  const char* operator*() const requires (not use_heap)
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

  char getc()
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
  using buf_type = std::conditional_t<use_heap, std::string, char[buf_size]>;
  using out_type = std::conditional_t<use_heap, std::string*, char*>;

  GttlFpType file;
  out_type out;
  bool is_end;
  buf_type default_buffer;
  size_t line_number;

  const char* input_string = nullptr;
  size_t string_length = 0;
  const char* current_ptr = nullptr;
};

#endif  // GTTL_LINE_GENERATOR_HPP
