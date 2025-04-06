#ifndef GTTL_LINE_GENERATOR_HPP
#define GTTL_LINE_GENERATOR_HPP

#include "utilities/gttl_file_open.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <iostream>
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
    if(length != nullptr) *length = 0;
    size_t len = 0;

    if constexpr(use_heap) out->clear();

    while(true)
    {
      if(file_buf_pos >= file_buf_end)
      {
        if(not refill_file_buffer())
        {
          is_end = true;
          return len > 0;
        }
      }

      char* next_newline =
        static_cast<char*>(
          std::memchr(file_buf + file_buf_pos,
                      '\n',
                      file_buf_end - file_buf_pos));
      if(next_newline != nullptr)
      {
        size_t line_len = next_newline - (file_buf + file_buf_pos);

        bool removed_cr = false;
        if(line_len > 0)
        {
          //Remove \r in case of DOS-style line-endings
          const char* const last_char = (file_buf + file_buf_pos + line_len - 1);
          if (*last_char == '\r')
          {
            line_len--;
            removed_cr = true;
          }
        }

        if constexpr(use_heap)
        {
          out->append(file_buf + file_buf_pos, line_len);
        }else
        {
          std::memcpy(out + len, file_buf + file_buf_pos, line_len);
          len += line_len;
          out[len] = '\0';
        }
        file_buf_pos += line_len + 1 + static_cast<size_t>(removed_cr);
        break;
      }
      // In this case there is no newline,
      // so we copy everything and continue reading
      size_t remaining = file_buf_end - file_buf_pos;
      if constexpr (use_heap)
      {
        out->append(file_buf + file_buf_pos, remaining);
      }else
      {
        std::memcpy(out + len, file_buf + file_buf_pos, remaining);
        len += remaining;
      }
      file_buf_pos = file_buf_end;
    }

    if(length != nullptr)
    {
      *length = len;
    }
    return true;
  }

  bool discard_line(size_t *length = nullptr)
  {
    char ch = EOF;
    size_t len = 0;

    if(length != nullptr)
    {
      *length = 0;
    }

    while(true)
    {
      if(file_buf_pos >= file_buf_end)
      {
        if(not refill_file_buffer())
        {
          is_end = true;
          return len > 0;
        }
      }

      ch = file_buf[file_buf_pos];
      if(ch == '\n')
      {
        file_buf_pos++;
        break;
      }

      file_buf_pos++;
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
    file_buf_pos = 0;
    file_buf_end = 0;
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

  void set_out_buffer(char* const _out) requires (not use_heap)
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
      if(current_ptr >= input_string + string_length)
      {
        is_end = true;
        return EOF;
      }
      return *current_ptr++;
    }

    if(file_buf_pos >= file_buf_end)
    {
      if(not refill_file_buffer())
      {
        is_end = true;
        return EOF;
      }
    }
    return file_buf[file_buf_pos++];
  }

  private:
  using buf_type = std::conditional_t<use_heap, std::string, char[buf_size]>;
  using out_type = std::conditional_t<use_heap, std::string*, char*>;

  char file_buf[buf_size]{};
  size_t file_buf_pos = 0;
  size_t file_buf_end = 0;

  bool refill_file_buffer()
  {
    file_buf_end = gttl_fp_type_read(file_buf, sizeof(char), buf_size, file);
    file_buf_pos = 0;
    return file_buf_end > 0;
  }

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
