#ifndef GTTL_LINE_GENERATOR_HPP
#define GTTL_LINE_GENERATOR_HPP

#include "utilities/gttl_file_open.hpp"
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <ios>
#include <span>
#include <string>
#include <vector>

template <const size_t buf_size = (size_t{1} << size_t{14}),
          const bool skip_empty_lines = false>
class GttlLineGenerator
{
  public:
  using CharBuffer = std::array<char, buf_size>;
  using BufferSpan = std::span<char>;

  explicit GttlLineGenerator (GttlFpType fp, bool _is_end = false)
    : file(fp)
    , out(&default_buffer)
    , is_end(_is_end)
    , line_number(1)
  {
    if(file == nullptr)
    {
      throw std::ios_base::failure(": cannot open file");
    }
  }

  explicit GttlLineGenerator(const char* file_name, bool _is_end = false)
    : file(gttl_fp_type_open(file_name, "rb"))
    , out(&default_buffer)
    , is_end(_is_end)
    , line_number(1)
  {
    if(file == nullptr)
    {
      throw std::ios_base::failure(": cannot open file");
    }
  }

  explicit GttlLineGenerator(const std::string& file_name, bool _is_end = false)
    : file(gttl_fp_type_open(file_name.c_str(), "rb"))
    , out(&default_buffer)
    , is_end(_is_end)
    , line_number(1)
  {
    if(file == nullptr)
    {
      throw std::ios_base::failure(": cannot open file");
    }
  }

  explicit GttlLineGenerator(const char* file_name,
                    std::string* _out,
                    bool _is_end = false)
    : file(gttl_fp_type_open(file_name, "rb"))
    , out(_out)
    , is_end(_is_end)
    , line_number(1)
  {
    if (file == nullptr)
    {
      throw std::ios_base::failure(": cannot open file");
    }
  }

  explicit GttlLineGenerator(const std::string& file_name,
                    std::string* _out,
                    bool _is_end = false)
    : file(gttl_fp_type_open(file_name.c_str(), "rb"))
    , out(_out)
    , is_end(_is_end)
    , line_number(1)
  {
    if (file == nullptr)
    {
      throw std::ios_base::failure(": cannot open file");
    }
  }

  explicit GttlLineGenerator(GttlFpType fp,
                             std::string* _out,
                             bool _is_end = false)
    : file(fp)
    , out(_out)
    , is_end(_is_end)
    , line_number(1)
  {
    if (file == nullptr)
    {
      throw std::ios_base::failure(": cannot open file");
    }
  }

  explicit GttlLineGenerator(const char* _input_string, size_t _string_length)
    : file(nullptr)
    , out(&default_buffer)
    , is_end(_string_length == 0)
    , line_number(1)
    , input_view(_input_string, _string_length)
    , current_ptr(input_view.data())
    {}

  explicit GttlLineGenerator(const std::vector<std::string>* _file_list,
                             std::string* _out = nullptr,
                             bool _is_end = false)
    : file(nullptr)
    , out(_out == nullptr ? &default_buffer : _out)
    , is_end(_file_list == nullptr or _file_list->empty() or _is_end)
    , line_number(1)
    , file_list(_file_list)
  {
    if (not is_end)
    {
      file = gttl_fp_type_open((*file_list)[file_index].c_str(), "rb");
      if (file == nullptr)
      {
        throw std::ios_base::failure(": cannot open file");
      }
    }
  }

  // Delete copy/move constructur & assignment operator
  GttlLineGenerator(const GttlLineGenerator&) = delete;
  GttlLineGenerator& operator=(const GttlLineGenerator&) = delete;
  GttlLineGenerator(GttlLineGenerator&&) = delete;
  GttlLineGenerator& operator=(GttlLineGenerator&&) = delete;

  ~GttlLineGenerator()
  {
    gttl_fp_type_close(file);
  }

private:
  bool read_from_mapped_string(size_t* length, bool append)
  {
    if (current_ptr >= input_view.data() + input_view.size())
    {
      is_end = true;
      return false;
    }

    const char* const next_newline = std::find(current_ptr,
                                       input_view.data() + input_view.size(),
                                       '\n');
    const size_t line_len = next_newline - current_ptr;
    const size_t copy_len = std::min(line_len, buf_size - 1);

    if (out != nullptr)
    {
      if (append)
      {
        const size_t old_size = out->size();
        out->resize(old_size + copy_len);
        std::memcpy(out->data() + old_size, current_ptr, copy_len);
      }
      else
      {
        out->resize(copy_len);
        std::memcpy(out->data(), current_ptr, copy_len);
      }
    }

    if (length != nullptr)
    {
      *length = copy_len;
    }

    current_ptr = (next_newline < input_view.data() + input_view.size())
                    ? next_newline + 1
                    : input_view.data() + input_view.size();
    return true;
  }

  bool read_from_file(size_t* length)
  {
    if (length != nullptr) *length = 0;
    size_t len = 0;

    while (true)
    {
      if (file_buf_pos >= file_buf_end)
      {
        if (not refill_file_buffer())
        {
          is_end = true;
          return len > 0;
        }
      }

      char* const next_newline = static_cast<char*>(std::memchr(
                                                      file_buf_span.data() +
                                                        file_buf_pos,
                                                      '\n',
                                                      file_buf_end -
                                                        file_buf_pos));
      if (next_newline != nullptr)
      {
        size_t line_len = next_newline - (file_buf_span.data() + file_buf_pos);

        bool removed_cr = false;
        if (line_len > 0 && file_buf_span[file_buf_pos + line_len - 1] == '\r')
        {
          line_len--;
          removed_cr = true;
        }

        out->append(file_buf_span.data() + file_buf_pos, line_len);
        len += line_len;
        file_buf_pos += line_len + 1 + static_cast<size_t>(removed_cr);
        break;
      }

      const size_t remaining = file_buf_end - file_buf_pos;
      out->append(file_buf_span.data() + file_buf_pos, remaining);
      len += remaining;
      file_buf_pos = file_buf_end;
    }

    if (length != nullptr)
    {
      *length = len;
    }
    return true;
  }

  bool discard_line(size_t* length = nullptr)
  {
    char ch = EOF;
    const size_t len = 0;

    if (length != nullptr)
    {
      *length = 0;
    }

    while (true)
    {
      if (file_buf_pos >= file_buf_end)
      {
        if (not refill_file_buffer())
        {
          is_end = true;
          return len > 0;
        }
      }

      ch = file_buf_span[file_buf_pos];
      if (ch == '\n')
      {
        file_buf_pos++;
        break;
      }

      file_buf_pos++;
      if (length != nullptr)
      {
        (*length)++;
      }
    }
    return (ch != EOF);
  }

public:
  bool advance(size_t* length = nullptr, bool append = false)
  {
    if (length != nullptr)
    {
      *length = 0;
    }
    if (is_end) return false;
    ++line_number;

    if (out != nullptr && !append) out->clear();

    size_t local_len = 0;
    bool ok = false;

    if (!input_view.empty())
    {
      ok = read_from_mapped_string(&local_len, append);
    }
    else if (out == nullptr)
    {
      ok = discard_line(&local_len);
    }
    else
    {
      ok = read_from_file(&local_len);
    }

    if (!ok) return false;

    if constexpr (skip_empty_lines)
    {
      if (local_len == 0 && !is_end && !line_partly_read)
      {
        return advance(length, append);
      }
      line_partly_read = false;
    }

    if (length != nullptr)
    {
      *length = local_len;
    }
    return true;
  }

  void reset()
  {
    is_end = false;
    line_number = 1;
    if (!input_view.empty())
    {
      current_ptr = input_view.data();
    }
    file_buf_pos = 0;
    file_buf_end = 0;
    if (file_list != nullptr && !file_list->empty())
    {
      gttl_fp_type_close(file);
      file_index = 0;
      file = gttl_fp_type_open((*file_list)[0].c_str(), "rb");
    }
    else
    {
      gttl_fp_type_rewind(file);
    }
  }

  [[nodiscard]] size_t line_number_get() const noexcept
  {
    return line_number;
  }

  void set_out_buffer(std::string* _out)
  {
    out = _out;
  }

  class Iterator
  {
  public:
    explicit Iterator(GttlLineGenerator* generator, bool end = false)
      : gen(generator), is_end(end)
    {
      if (!end) ++(*this);
    }

    const std::string& operator*() const
    {
      return *(gen->out);
    }

    Iterator& operator++()
    {
      is_end = !gen->advance();
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

  std::string operator*() const
  {
    return *out;
  }

  Iterator begin() { return Iterator(this, false); }
  Iterator end() { return Iterator(this, true); }

  char getc()
  {
    if (!input_view.empty())
    {
      if (current_ptr >= input_view.data() + input_view.size())
      {
        is_end = true;
        return EOF;
      }
      if constexpr (skip_empty_lines) line_partly_read = true;
      return *current_ptr++;
    }

    if (file_buf_pos >= file_buf_end)
    {
      if (!refill_file_buffer())
      {
        is_end = true;
        return EOF;
      }
    }
    if constexpr (skip_empty_lines) line_partly_read = true;
    return file_buf_span[file_buf_pos++];
  }

private:
  CharBuffer file_buf{};
  BufferSpan file_buf_span{file_buf};
  size_t file_buf_pos = 0;
  size_t file_buf_end = 0;

  bool refill_file_buffer()
  {
    while (true)
    {
      if (file == nullptr)
      {
        is_end = true;
        return false;
      }
      file_buf_end = gttl_fp_type_read(file_buf.data(),
                                       sizeof(char),
                                       buf_size,
                                       file);
      file_buf_pos = 0;
      if (file_buf_end > 0)
      {
        return true;
      }
      if (file_list == nullptr)
      {
        is_end = true;
        return false;
      }
      gttl_fp_type_close(file);
      file = nullptr;

      if (++file_index < file_list->size())
      {
        file = gttl_fp_type_open((*file_list)[file_index].c_str(), "rb");
        if (file == nullptr)
        {
          throw std::ios_base::failure(": cannot open file");
        }
      }
    }
  }

  GttlFpType file;
  std::string* out;
  bool is_end;
  std::string default_buffer;
  size_t line_number;

  std::span<const char> input_view;
  const char* current_ptr = nullptr;

  const std::vector<std::string>* file_list = nullptr;
  size_t file_index = 0;

  bool line_partly_read = false;
};

#endif // GTTL_LINE_GENERATOR_HPP
