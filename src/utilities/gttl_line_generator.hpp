#ifndef GTTL_LINE_GENERATOR_HPP
#define GTTL_LINE_GENERATOR_HPP

#include "utilities/gttl_file_open.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <ios>
#include <span>
#include <string>
#include <utility>
#include <vector>

template <size_t buf_size = (size_t{1} << 14),
          bool skip_empty_lines = false>
class GttlLineGenerator
{
  using CharBuffer = std::array<char, buf_size>;
  using BufferSpan = std::span<char>;

  // fixed-size buffer for reading a chunk of data from the file
  CharBuffer file_buf{};
  // A std::span view into the buffer, to provide abstraction over the
  // otherwise complex and error-prone pointer arithmetic
  BufferSpan file_buf_span{file_buf};
  // Current read position within the file_buf buffer
  size_t file_buf_pos = 0;
  // index marking the end of valid content in file_buf, ie.
  // the size of the buffer's actual content plus one
  size_t file_buf_end = 0;

  // Handle to the currently open file from which we are reading
  GttlFpType file;
  // output-pointer towards the string in which the current line should be
  // stored.
  std::string* line_ptr;
  // default output-buffer for when we do not want to read directly into a
  // different object's memory
  std::string default_buffer;
  // flag to indicate whether we are fully done reading
  bool all_files_exhausted;
  // current line counter. Incremented whenever we read a line
  size_t line_number;

  // This is a view into a memory-mapped input string, to simplify
  // pointer-arithmetic (as done above for files) when we read from
  // such a string.
  std::span<const char> input_view;
  // And the current read-position within it, analogous to file_buf_pos
  const char* current_ptr = nullptr;

  // A vector of multiples files from which we may want to read.
  // Allocated only when the appropriate constructor is called
  const std::vector<std::string>* file_list = nullptr;
  // The index of the currently read file.
  // Files will be read consecutively, as though concatenated into a single file
  size_t file_index = 0;

  /* this is only used skip_empty_lines is true */
  bool line_partly_read = false;

  bool refill_file_buffer(void)
  {
    while (true)
    {
      if (file == nullptr)
      {
        all_files_exhausted = true;
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
        all_files_exhausted = true;
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

  bool read_from_mapped_string(size_t* length_ptr, bool append)
  {
    assert(length_ptr != nullptr and *length_ptr == 0);
    if (current_ptr >= input_view.data() + input_view.size())
    {
      all_files_exhausted = true;
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
          all_files_exhausted = true;
          return len > 0;
        }
      }

      const char* const next_newline
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
          all_files_exhausted = true;
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
  explicit GttlLineGenerator (GttlFpType fp, bool _exhausted = false)
    : file(fp)
    , line_ptr(&default_buffer)
    , all_files_exhausted(_exhausted)
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
    , all_files_exhausted(_exhausted)
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
    , all_files_exhausted(_exhausted)
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
    , all_files_exhausted(_exhausted)
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
    , all_files_exhausted(_exhausted)
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
    , all_files_exhausted(_exhausted)
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
    , all_files_exhausted(_string_length == 0 or _exhausted)
    , line_number(1)
    , input_view(_input_string, _string_length)
    , current_ptr(input_view.data())
  { }

  explicit GttlLineGenerator(const std::vector<std::string>* _file_list,
                             std::string* _line_ptr = nullptr,
                             bool _exhausted = false)
    : file(nullptr)
    , line_ptr(_line_ptr == nullptr ? &default_buffer : _line_ptr)
    , all_files_exhausted(_file_list == nullptr or _file_list->empty()
                                                or _exhausted)
    , line_number(1)
    , file_list(_file_list)
  {
    if (not all_files_exhausted)
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

  std::pair<bool, size_t> advance(bool append = false)
  {
    if (all_files_exhausted)
    {
      return std::make_pair(false, 0);
    }
    ++line_number;

    if (line_ptr != nullptr and not append)
    {
      line_ptr->clear();
    }

    size_t local_len = 0;
    bool okay;

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
      return std::make_pair(false, local_len);
    }

    if constexpr (skip_empty_lines)
    {
      if (local_len == 0 and not all_files_exhausted and not line_partly_read)
      {
        return advance(append);
      }
      line_partly_read = false;
    }

    return std::make_pair(true, local_len);
  }

  char getc(void)
  {
    if (not input_view.empty())
    {
      if (current_ptr >= input_view.data() + input_view.size())
      {
        all_files_exhausted = true;
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
        all_files_exhausted = true;
        return EOF;
      }
    }
    if constexpr (skip_empty_lines)
    {
      line_partly_read = true;
    }
    return file_buf_span[file_buf_pos++];
  }

  void reset(void)
  {
    all_files_exhausted = false;
    line_number = 1;
    if (not input_view.empty())
    {
      // current_ptr is only relevant when we read from a memory-mapped string.
      // Hence why we only reset it when there exists an input_view
      // When reading from a file, we instead use file_buf_pos/file_buf_end
      current_ptr = input_view.data();
    }
    file_buf_pos = 0;
    file_buf_end = 0;
    if (file_list != nullptr and not file_list->empty())
    {
      gttl_fp_type_close(file);
      file_index = 0;
      file = gttl_fp_type_open((*file_list)[0].c_str(), "rb");
    } else
    {
      gttl_fp_type_rewind(file);
    }
    // line_ptr may point to an externally provided memory buffer, hence
    // why it isn't reset.
    // Setting it to an empty string may overwrite data that would be
    // accessed elsewhere. We are primarily providing data to be processed in
    // other places, and copying or moving ownership here would be an additional
    // cost. Hence the set_line_buffer() option.
    //
    // line_partly_read is simply irrelevant after a reset,
    // since we begin reading anew
    // regardless.
    // It will be overwritten on the first call to advance()
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

  [[nodiscard]] const std::string& data_get(void) const
  {
    return *line_ptr;
  }

  class Iterator
  {
    private:
    GttlLineGenerator* generator;
    bool iter_exhausted;
    public:
    explicit Iterator(GttlLineGenerator* _generator, bool _exhausted)
      : generator(_generator)
      , iter_exhausted(_exhausted)
    {
      if (not iter_exhausted)
      {
        ++(*this);
      }
    }

    auto operator * (void) const
    {
      return generator->data_get();
    }

    const Iterator& operator ++ (void)
    {
      assert(not iter_exhausted);
      iter_exhausted = not std::get<0>(generator->advance());
      return *this;
    }

    bool operator == (const Iterator& other) const
    {
      return iter_exhausted == other.iter_exhausted and
             generator == other.generator;
    }

    bool operator != (const Iterator& other) const
    {
      return not (*this == other);
    }
  };

  Iterator begin(void)
  {
    return Iterator(this, false);
  }
  Iterator end(void)
  {
    return Iterator(this, true);
  }
};
#endif // GTTL_LINE_GENERATOR_HPP
