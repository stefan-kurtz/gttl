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
    //TODO: Test this case. In all generators.
    if(out == nullptr)
    {
      if constexpr(use_heap)
      {
        std::string to_discard;
        lg.set_out_buffer(&to_discard);
      }else
      {
        char to_discard[buf_size];
        lg.set_out_buffer(to_discard);
      }
      lg.advance();
      lg.advance();
      is_first_entry = false;
      while(lg.getc() != '>')
      {
        if(not lg.advance()) return false;
      }
      return true;
    }

    if constexpr (use_heap)
    {
      out->header.clear();
      out->sequence.clear();
      if (is_first_entry)
      {
        lg.getc();
        is_first_entry = false;
      }
      lg.set_out_buffer(&out->header);
      if(not lg.advance()) return false;
    }else
    {
      if(is_first_entry)
      {
        lg.getc();
        is_first_entry = false;
      }
      lg.set_out_buffer(out->header);
      if(not lg.advance()) return false;
    }
    size_t offset = 0;
    size_t line_length = 0;
    char next = lg.getc();
    // (char) -1 is EOF
    while(next != '>' and next != static_cast<char>(-1))
    {
      if constexpr(use_heap)
      {
        lg.set_out_buffer(&out->sequence);
        lg.advance(&line_length);
        out->sequence.insert(0, 1, next);
      }else
      {
        out->sequence[offset++] = next;
        lg.set_out_buffer((out->sequence) + ((offset) * sizeof(char)));
        lg.advance(&line_length);
        offset += line_length;
      }
      next = lg.getc();
    }
    return true;
  }

  void reset(void)
  {
    lg.reset();
    is_first_entry = true;
    is_end = false;
  }

  void line_number(void) const noexcept
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
