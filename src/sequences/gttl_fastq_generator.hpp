#ifndef GTTL_FASTQ_GENERATOR_HPP
#define GTTL_FASTQ_GENERATOR_HPP

#include <cstddef>
#include <string_view>
#include "utilities/gttl_file_open.hpp"
#include "utilities/gttl_line_generator.hpp"

template <const size_t buf_size = (1 << 14)>
struct GttlFastQEntry
{
  std::string_view header_get() const noexcept
  {
    return std::string_view(header);
  }
  std::string_view sequence_get() const noexcept
  {
    return std::string_view(sequence);
  }
  std::string_view quality_get() const noexcept
  {
    return std::string_view(quality);
  }
  char header[buf_size];
  char sequence[buf_size];
  char quality[buf_size];
};

template <const size_t buf_size = (1 << 14)>
class GttlFastQGenerator
{
  public:
  GttlFastQGenerator(const char* file_name,
                     bool _is_end = false)
    : out(&default_buffer)
    , is_end(_is_end)
    , lg(gttl_fp_type_open(file_name, "rb"), out->header, _is_end) {}

  GttlFastQGenerator(GttlFpType fp,
                     bool _is_end = false)
    : out(&default_buffer)
    , is_end(_is_end)
    , lg(fp, out->header, _is_end) {}


  GttlFastQGenerator(const char* file_name,
                     GttlFastQEntry<buf_size> *_out,
                     bool _is_end = false)
    : lg(gttl_fp_type_open(file_name, "rb"), _out->header, _is_end)
    , out(_out)
    , is_end(_is_end) {}

  GttlFastQGenerator(GttlFpType fp,
                     GttlFastQEntry<buf_size> *_out,
                     bool _is_end = false)
    : lg(fp, _out->header, _is_end)
    , out(_out)
    , is_end(_is_end) {}

  GttlFastQGenerator(const char* _input_string,
                     size_t _string_length,
                     bool _is_end = false)
    : out(&default_buffer)
    , is_end(_is_end)
    , lg(_input_string, _string_length)
    {}

  bool advance(void)
  {
    if(is_end) return false;
    if(out == nullptr) return true;
    lg.set_out_buffer(out->header);
    lg.advance();
    lg.set_out_buffer(out->sequence);
    lg.advance();
    lg.set_out_buffer(nullptr);
    lg.advance();
    lg.set_out_buffer(out->quality);
    return lg.advance();
  }

  void reset(void)
  {
    lg.reset();
  }

  size_t line_number(void) const noexcept
  {
    return lg.line_number_get();
  }

  class Iterator
  {
    public:
    Iterator(GttlFastQGenerator* generator, bool end = false)
      : gen(generator), is_end(end)
    {
      if(!end) ++(*this);
    }

    const GttlFastQEntry<buf_size>* operator*() const
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
    GttlFastQGenerator* gen;
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
  GttlFastQEntry<buf_size> default_buffer;
  GttlFastQEntry<buf_size>* out;
  bool is_end;
  GttlLineGenerator<buf_size> lg;
};

#endif  // GTTL_FASTQ_GENERATOR_HPP
