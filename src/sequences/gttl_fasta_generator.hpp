#ifndef GTTL_FASTA_GENERATOR_HPP
#define GTTL_FASTA_GENERATOR_HPP

#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include "utilities/gttl_file_open.hpp"
#include "utilities/gttl_line_generator.hpp"

template <const size_t buf_size = (1 << 14), const bool use_heap = false>
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

template <const size_t buf_size = (1 << 14), const bool use_heap = false>
class GttlFastAGenerator
{
  public:
  explicit GttlFastAGenerator(const char* file_name,
                              bool _is_end = false) requires (!use_heap)
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
      //TODO: Test this! Both here and in the FastQGenerator.
      // This section might still be affected by the return-value problem.
      if constexpr(use_heap)
      {
        std::string to_discard;
        lg.set_out_buffer(&to_discard);
        lg.advance();
        lg.advance();
        bool ret = true;
        while(!to_discard.starts_with('>'))
        {
          ret = lg.advance();
        }
        return ret;
      }else
      {
        char to_discard[buf_size];
        lg.set_out_buffer(to_discard);
        lg.advance();
        lg.advance();
        bool ret = true;
        while(to_discard[0] != '>')
        {
          ret = lg.advance();
        }
        return ret;
      }
    }
    char next = lg.getc();
    // (char) -1 is EOF
    if(next == -1) return false;

    if constexpr (use_heap)
    {
      out->header.clear();
      out->sequence.clear();
      if (next == '>')
      {
        lg.set_out_buffer(&out->header);
        lg.advance();
        out->header = '>' + out->header;
      }else
      {
        lg.set_out_buffer(&out->header);
        lg.advance();
        out->header = ">";
        out->header += next + out->header;
      }
    }else
    {
      out->header[0] = '>';
      if (next == '>')
      {
        lg.set_out_buffer(out->header + sizeof(char));
      } else
      {
        (out->header)[1] = next;
        lg.set_out_buffer(out->header + (2 * sizeof(char)));
      }
      lg.advance();
    }
    size_t offset = 0;
    size_t line_length = 0;
    next = lg.getc();
    while(next != '>')
    {
      // (char) -1 is EOF
      if(next == -1)
        break;
      if constexpr(use_heap)
      {
        //TODO: Is there a better way to do this?
        //      frequently copying a buffer like this is slow.
        std::string next_seqline;
        lg.set_out_buffer(&next_seqline);
        lg.advance(&line_length);
        out->sequence += next + next_seqline;
        offset += line_length + 1;
      }else
      {
        out->sequence[offset] = next;
        offset++;
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
      if(!end) ++(*this);
    }

    const GttlFastAEntry<buf_size, use_heap>* operator*() const
    {
      return gen->out;
    }

    const Iterator& operator++()
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
  GttlLineGenerator<buf_size, use_heap> lg;
};
#endif // GTTL_FASTA_GENERATOR_HPP
