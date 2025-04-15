#ifndef GTTL_LINE_GENERATOR_HPP
#define GTTL_LINE_GENERATOR_HPP

#include "utilities/gttl_file_open.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <string>

template <const size_t buf_size = (size_t{1} << size_t{14})>
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
                    bool _is_end = false)
    : file(gttl_fp_type_open(file_name, "rb"))
    , out(_out)
    , is_end(_is_end)
    , line_number(0)
    {}

  explicit GttlLineGenerator(GttlFpType fp,
                             std::string* _out,
                             bool _is_end = false)
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
      out->resize(copy_len);
      std::memcpy(out->data(), current_ptr, copy_len);
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
          const char* const last_char =
            (file_buf + file_buf_pos + line_len - 1);
          if (*last_char == '\r')
          {
            line_len--;
            removed_cr = true;
          }
        }

        out->append(file_buf + file_buf_pos, line_len);
        file_buf_pos += line_len + 1 + static_cast<size_t>(removed_cr);
        break;
      }
      // In this case there is no newline,
      // so we copy everything and continue reading
      size_t remaining = file_buf_end - file_buf_pos;
      out->append(file_buf + file_buf_pos, remaining);
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
  bool advance(size_t *length = nullptr, bool append = false)
  {
    if(length != nullptr)
    {
      *length = 0;
    }
    if(is_end) return false;
    ++line_number;

    if(out != nullptr and not append) out->clear();

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

  void set_out_buffer(std::string* _out)
  {
    out = _out;
  }

  void set_default_out_buffer()
  {
    out = &default_buffer;
  }

  class Iterator
  {
    public:
    explicit Iterator(GttlLineGenerator* generator, bool end = false)
      : gen(generator), is_end(end)
    {
      if(not end) ++(*this);
    }

    const std::string& operator*() const
    {
      return *(gen->out);
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

  std::string operator*() const
  {
    return *out;
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
  std::string* out;
  bool is_end;
  std::string default_buffer;
  size_t line_number;

  const char* input_string = nullptr;
  size_t string_length = 0;
  const char* current_ptr = nullptr;
};

#endif  // GTTL_LINE_GENERATOR_HPP
