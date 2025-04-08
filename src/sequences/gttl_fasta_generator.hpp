#ifndef GTTL_FASTA_GENERATOR_HPP
#define GTTL_FASTA_GENERATOR_HPP

#include "utilities/gttl_file_open.hpp"
#include "utilities/gttl_line_generator.hpp"
#include <cstddef>
#include <string>
#include <string_view>
#include <type_traits>

template <const size_t buf_size = (1U << 14U), const bool use_heap = false>
struct GttlFastAEntry
{
  [[nodiscard]] std::string_view header_get() const noexcept
  {
    return std::string_view(header);
  }
  [[nodiscard]] std::string_view sequence_get() const noexcept
  {
    return std::string_view(sequence);
  }

  using buf_type = std::conditional_t<use_heap, std::string, char[buf_size]>;
  buf_type header;
  buf_type sequence;

  GttlFastAEntry()
  {
    if constexpr (use_heap)
    {
      header.reserve(buf_size);
      sequence.reserve(buf_size);
    }
  }
};

template <const size_t buf_size = (1U << 14U), const bool use_heap = false>
class GttlFastAGenerator
{
  public:
  explicit GttlFastAGenerator(const char* file_name,
                              bool _is_end = false) requires (not use_heap)
    : out(&default_buffer)
    , is_end(_is_end)
    , lg(gttl_fp_type_open(file_name, "rb"), out->header, _is_end)
  {}

  explicit GttlFastAGenerator(const char* file_name,
                              bool _is_end = false) requires use_heap
    : out(&default_buffer)
    , is_end(_is_end)
    , lg(gttl_fp_type_open(file_name, "rb"), &(out->header), _is_end)
  {}

  explicit GttlFastAGenerator(GttlFpType fp,
                              bool _is_end = false)
    : out(&default_buffer)
    , is_end(_is_end)
    , lg(fp, out->header, _is_end)
  {}

  explicit GttlFastAGenerator(const char* file_name,
                              GttlFastAEntry<buf_size, use_heap> *_out,
                              bool _is_end = false)
    : out(_out)
    , is_end(_is_end)
    , lg(gttl_fp_type_open(file_name, "rb"), _out->header, _is_end)
  {}

  explicit GttlFastAGenerator(GttlFpType fp,
                              GttlFastAEntry<buf_size, use_heap> *_out,
                              bool _is_end = false)
    : out(_out)
    , is_end(_is_end)
    , lg(fp, _out->header, _is_end)
  {}

  explicit GttlFastAGenerator(const char* _input_string,
                              size_t _string_length,
                              bool _is_end = false)
    : out(&default_buffer)
    , is_end(_is_end)
    , lg(_input_string, _string_length)
  {}

  bool advance()
  {
    if(is_end) return false;

    if(out == nullptr)
    {
      lg.set_out_buffer(nullptr);

      if(not lg.advance()) return false;

      while(true)
      {
        int ch = lg.getc();
        if(ch == EOF or ch == '>') break;
        lg.advance();
      }

      is_end = (lg.line_number_get() == 0);
      return not is_end;
    }

    if(is_first_entry)
    {
      int ch = lg.getc();
      if(ch != '>') return false;
      is_first_entry = false;
    }

    if constexpr(use_heap)
    {
      out->header.clear();
      out->sequence.clear();
      lg.set_out_buffer(&out->header);
    }else
    {
      out->header[0] = '\0';
      out->sequence[0] = '\0';
      lg.set_out_buffer(out->header);
    }

    if(not lg.advance())
    {
      is_end = true;
      return false;
    }

    if constexpr(not use_heap)
    {
      size_t offset = 0;
      while(true)
      {
        int ch = lg.getc();
        if(ch == EOF or ch == '>') break;

        out->sequence[offset++] = static_cast<char>(ch);
        lg.set_out_buffer(out->sequence + offset);
        size_t len = 0;
        if(not lg.advance(&len)) break;
        offset += len;
      }
      out->sequence[offset] = '\0';
    }else
    {
      lg.set_out_buffer(&out->sequence);
      while(true)
      {
        int ch = lg.getc();
        if(ch == EOF or ch == '>') break;

        out->sequence.push_back(static_cast<char>(ch));
        if(not lg.advance(nullptr, true)) break;
      }
    }
    return true;
  }

  void reset(void)
  {
    lg.reset();
    is_first_entry = true;
    is_end = false;
  }

  [[nodiscard]] size_t line_number(void) const noexcept
  {
    return lg.line_number_get();
  }

  class Iterator
  {
    public:
    explicit Iterator(GttlFastAGenerator* generator, bool end = false)
      : gen(generator), is_end(end)
    {
      if(not end) ++(*this);
    }

    const GttlFastAEntry<buf_size, use_heap>* operator*() const
    {
      return gen->out;
    }

    const Iterator& operator++()
    {
      if(not gen->advance()) is_end = true;
      return *this;
    }

    bool operator==(const Iterator& other) const
    {
      return (is_end == other.is_end) and (gen == other.gen);
    }

    bool operator!=(const Iterator& other) const
    {
      return not (*this == other);
    }

    private:
    GttlFastAGenerator* gen;
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
  GttlFastAEntry<buf_size, use_heap> default_buffer;
  GttlFastAEntry<buf_size, use_heap>* out;
  bool is_end;
  bool is_first_entry = true;
  GttlLineGenerator<buf_size, use_heap> lg;
};
#endif // GTTL_FASTA_GENERATOR_HPP
