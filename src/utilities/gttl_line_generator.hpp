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

template <size_t buf_size = (size_t{1} << size_t{14}),
          bool skip_empty_lines = false>
class GttlLineGenerator
{
  private:
  using CharBuffer = std::array<char, buf_size>;
  using BufferSpan = std::span<char>;
  CharBuffer file_buf{};
  BufferSpan file_buf_span{file_buf};
  size_t file_buf_pos = 0;
  size_t file_buf_end = 0;

  GttlFpType file;
  std::string* line_ptr;
  bool exhausted;
  std::string default_buffer;
  size_t line_number;

  std::span<const char> input_view;
  const char* current_ptr = nullptr;

  const std::vector<std::string>* file_list = nullptr;
  size_t file_index = 0;

  bool line_partly_read = false;

  bool refill_file_buffer(void)
  {
    while (true)
    {
      if (file == nullptr)
      {
        exhausted = true;
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
        exhausted = true;
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

  public:
  explicit GttlLineGenerator (GttlFpType fp, bool _exhausted = false)
    : file(fp)
    , line_ptr(&default_buffer)
    , exhausted(_exhausted)
    , line_number(1)
  {
    if (file == nullptr)
    {
      throw std::ios_base::failure(": cannot open file");
    }
  }

  explicit GttlLineGenerator(const char* file_name, bool _exhausted = false)
    : file(gttl_fp_type_open(file_name, "rb"))
    , line_ptr(&default_buffer)
    , exhausted(_exhausted)
    , line_number(1)
  {
    if (file == nullptr)
    {
      throw std::ios_base::failure(": cannot open file");
    }
  }

  explicit GttlLineGenerator(const std::string& file_name,
                             bool _exhausted = false)
    : file(gttl_fp_type_open(file_name.c_str(), "rb"))
    , line_ptr(&default_buffer)
    , exhausted(_exhausted)
    , line_number(1)
  {
    if (file == nullptr)
    {
      throw std::ios_base::failure(": cannot open file");
    }
  }

  explicit GttlLineGenerator(const char* file_name,
                             std::string* _line_ptr, bool _exhausted = false)
    : file(gttl_fp_type_open(file_name, "rb"))
    , line_ptr(_line_ptr)
    , exhausted(_exhausted)
    , line_number(1)
  {
    if (file == nullptr)
    {
      throw std::ios_base::failure(": cannot open file");
    }
  }

  explicit GttlLineGenerator(const std::string& file_name,
                             std::string* _line_ptr, bool _exhausted = false)
    : file(gttl_fp_type_open(file_name.c_str(), "rb"))
    , line_ptr(_line_ptr)
    , exhausted(_exhausted)
    , line_number(1)
  {
    if (file == nullptr)
    {
      throw std::ios_base::failure(": cannot open file");
    }
  }

  explicit GttlLineGenerator(GttlFpType fp,
                             std::string* _line_ptr, bool _exhausted = false)
    : file(fp)
    , line_ptr(_line_ptr)
    , exhausted(_exhausted)
    , line_number(1)
  {
    if (file == nullptr)
    {
      throw std::ios_base::failure(": cannot open file");
    }
  }

  explicit GttlLineGenerator(const char* _input_string, size_t _string_length,
                             bool _exhausted = false)
    : file(nullptr)
    , line_ptr(&default_buffer)
    , exhausted(_string_length == 0 or _exhausted)
    , line_number(1)
    , input_view(_input_string, _string_length)
    , current_ptr(input_view.data())
  { }

  explicit GttlLineGenerator(const std::vector<std::string>* _file_list,
                             std::string* _line_ptr = nullptr,
                             bool _exhausted = false)
    : file(nullptr)
    , line_ptr(_line_ptr == nullptr ? &default_buffer : _line_ptr)
    , exhausted(_file_list == nullptr or _file_list->empty() or _exhausted)
    , line_number(1)
    , file_list(_file_list)
  {
    if (not exhausted)
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

  ~GttlLineGenerator(void)
  {
    gttl_fp_type_close(file);
  }

  private:
  bool read_from_mapped_string(size_t* length_ptr, bool append)
  {
    assert(length_ptr != nullptr and *length_ptr == 0);
    if (current_ptr >= input_view.data() + input_view.size())
    {
      exhausted = true;
      return false;
    }

    const char* const next_newline
      = std::find(current_ptr,input_view.data() + input_view.size(),'\n');
    const size_t line_len = static_cast<size_t>(next_newline - current_ptr);
    const size_t copy_len = std::min(line_len, buf_size - 1);

    if (line_ptr != nullptr)
    {
      if (append)
      {
        const size_t old_size = line_ptr->size();
        line_ptr->resize(old_size + copy_len);
        std::memcpy(line_ptr->data() + old_size, current_ptr, copy_len);
      } else
      {
        line_ptr->resize(copy_len);
        std::memcpy(line_ptr->data(), current_ptr, copy_len);
      }
    }

    *length_ptr = copy_len;

    current_ptr = (next_newline < input_view.data() + input_view.size())
                    ? next_newline + 1
                    : input_view.data() + input_view.size();
    return true;
  }

  bool read_from_file(size_t* length_ptr)
  {
    assert(length_ptr != nullptr and *length_ptr == 0);
    size_t len = 0;

    while (true)
    {
      if (file_buf_pos >= file_buf_end)
      {
        if (not refill_file_buffer())
        {
          exhausted = true;
          return len > 0;
        }
      }

      char* const next_newline
        = static_cast<char*>(std::memchr(file_buf_span.data() + file_buf_pos,
                                         '\n',
                                         file_buf_end - file_buf_pos));
      if (next_newline != nullptr)
      {
        size_t line_len = next_newline - (file_buf_span.data() + file_buf_pos);

        bool removed_cr = false;
        if (line_len > 0 and file_buf_span[file_buf_pos + line_len - 1] == '\r')
        {
          line_len--;
          removed_cr = true;
        }

        line_ptr->append(file_buf_span.data() + file_buf_pos, line_len);
        len += line_len;
        file_buf_pos += line_len + 1 + static_cast<size_t>(removed_cr);
        break;
      }

      const size_t remaining = file_buf_end - file_buf_pos;
      line_ptr->append(file_buf_span.data() + file_buf_pos, remaining);
      len += remaining;
      file_buf_pos = file_buf_end;
    }

    *length_ptr = len;
    return true;
  }

  bool discard_line(size_t* length_ptr)
  {
    char ch = EOF;

    assert(length_ptr != nullptr and *length_ptr == 0);
    while (true)
    {
      if (file_buf_pos >= file_buf_end)
      {
        if (not refill_file_buffer())
        {
          exhausted = true;
          return false;
        }
      }

      ch = file_buf_span[file_buf_pos];
      if (ch == '\n')
      {
        file_buf_pos++;
        break;
      }

      file_buf_pos++;
      (*length_ptr)++;
    }
    return (ch != EOF);
  }

public:
  bool advance(size_t* length_ptr = nullptr, bool append = false)
  {
    if (length_ptr != nullptr)
    {
      *length_ptr = 0;
    }
    if (exhausted)
    {
      return false;
    }
    ++line_number;

    if (line_ptr != nullptr and not append)
    {
      line_ptr->clear();
    }

    size_t local_len = 0;
    bool okay;

    /* Do not use pointer to local_len, but return the length value
       with the boolean value as pair. */
    if (not input_view.empty())
    {
      okay = read_from_mapped_string(&local_len, append);
    } else
    {
      if (line_ptr == nullptr)
      {
        okay = discard_line(&local_len);
      } else
      {
        okay = read_from_file(&local_len);
      }
    }
    if (not okay)
    {
      return false;
    }

    if constexpr (skip_empty_lines)
    {
      if (local_len == 0 and not exhausted and not line_partly_read)
      {
        return advance(length_ptr, append);
      }
      line_partly_read = false;
    }

    if (length_ptr != nullptr)
    {
      *length_ptr = local_len;
    }
    return true;
  }

  void reset(void)
  {
    exhausted = false;
    line_number = 1;
    if (not input_view.empty())
    {
      current_ptr = input_view.data();
    }
    file_buf_pos = 0;
    file_buf_end = 0;
    if (file_list != nullptr and not file_list->empty())
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

  void set_line_buffer(std::string* _line_ptr)
  {
    line_ptr = _line_ptr;
  }

  /*std::string operator*(void) const
  {
    return *line_ptr;
  }*/

  private:
  class Iterator
  {
    private:
    GttlLineGenerator* generator;
    bool exhausted;
    public:
    explicit Iterator(GttlLineGenerator* _generator, bool _exhausted)
      : generator(_generator)
      , exhausted(_exhausted)
    {
      if (not exhausted)
      {
        ++(*this);
      }
    }

    const std::string& operator*(void) const
    {
      return *(generator->line_ptr);
    }

    Iterator& operator++()
    {
      assert(not exhausted);
      exhausted = not generator->advance();
      return *this;
    }

    bool operator == (const Iterator& other) const
    {
      return exhausted == other.exhausted and
             generator == other.generator;
    }

    bool operator != (const Iterator& other) const
    {
      return not (*this == other);
    }
  };
  public:

  Iterator begin(void)
  {
    return Iterator(this, false);
  }
  Iterator end(void)
  {
    return Iterator(this, true);
  }

  char getc(void)
  {
    if (not input_view.empty())
    {
      if (current_ptr >= input_view.data() + input_view.size())
      {
        exhausted = true;
        return EOF;
      }
      if constexpr (skip_empty_lines)
      {
        line_partly_read = true;
      }
      return *current_ptr++;
    }

    if (file_buf_pos >= file_buf_end)
    {
      if (not refill_file_buffer())
      {
        exhausted = true;
        return EOF;
      }
    }
    if constexpr (skip_empty_lines)
    {
      line_partly_read = true;
    }
    return file_buf_span[file_buf_pos++];
  }
};
#endif // GTTL_LINE_GENERATOR_HPP
