#ifndef GTTL_FASTQ_GENERATOR_HPP
#define GTTL_FASTQ_GENERATOR_HPP

#include "utilities/gttl_file_open.hpp"
#include "utilities/gttl_line_generator.hpp"
#include <cstddef>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <string_view>

template <const size_t buf_size = (size_t{1} << size_t{14})>
struct GttlFastQEntry
{
  [[nodiscard]] std::string_view header_get() const noexcept
  {
    return std::string_view(header);
  }
  [[nodiscard]] std::string_view sequence_get() const noexcept
  {
    return std::string_view(sequence);
  }
  [[nodiscard]] std::string_view quality_get() const noexcept
  {
    return std::string_view(quality);
  }

  std::string header;
  std::string sequence;
  std::string quality;

  GttlFastQEntry()
  {
    header.reserve(buf_size);
    sequence.reserve(buf_size);
    quality.reserve(buf_size);
  }
};

template <const size_t buf_size = (size_t{1} << size_t{14})>
class GttlFastQGenerator
{
  public:
  static constexpr const bool is_fastq_generator = true;
  explicit GttlFastQGenerator(const char* file_name,
                              bool _is_end = false)
    : out(&default_buffer)
    , is_end(_is_end)
    , lg(gttl_fp_type_open(file_name, "rb"), &(out->header), _is_end) {}

  explicit GttlFastQGenerator(GttlFpType fp,
                              bool _is_end = false)
    : out(&default_buffer)
    , is_end(_is_end)
    , lg(fp, &out->header, _is_end) {}


  explicit GttlFastQGenerator(const char* file_name,
                              GttlFastQEntry<buf_size> *_out,
                              bool _is_end = false)
    : lg(gttl_fp_type_open(file_name, "rb"), _out->header, _is_end)
    , out(_out)
    , is_end(_is_end) {}

  explicit GttlFastQGenerator(GttlFpType fp,
                              GttlFastQEntry<buf_size> *_out,
                              bool _is_end = false)
    : lg(fp, _out->header, _is_end)
    , out(_out)
    , is_end(_is_end) {}

  explicit GttlFastQGenerator(const char* _input_string,
                              size_t _string_length,
                              bool _is_end = false)
    : out(&default_buffer)
    , is_end(_is_end)
    , lg(_input_string, _string_length)
    {}

  explicit GttlFastQGenerator(const std::vector<std::string>* _file_list,
                              GttlFastQEntry<buf_size> *_out = nullptr,
                              bool _is_end = false)
    : out(_out == nullptr ? &default_buffer : _out)
    , is_end(_is_end)
    , lg(_file_list, &out->header, _is_end)
    {}


  bool advance()
  {
    if(is_end) return false;
    if(out == nullptr)
    {
      lg.set_out_buffer(nullptr);
      for(size_t line = 0; line < 4; line++)
      {
        lg.advance();
      }
      return true;
    }
    lg.set_out_buffer(&out->header);
    const char c = lg.getc();
    if(c == EOF)
    {
      is_end = true;
      return false;
    }
    if(c != '@')
    {
      throw std::runtime_error(", line " + std::to_string(lg.line_number_get())
                               + ": corrupted sequence");
    }
    lg.advance();
    lg.set_out_buffer(&out->sequence);
    lg.advance();
    lg.set_out_buffer(nullptr);
    lg.advance();
    lg.set_out_buffer(&out->quality);
    return lg.advance();
  }

  void reset()
  {
    lg.reset();
    is_end = false;
  }

  [[nodiscard]] size_t line_number() const noexcept
  {
    return lg.line_number_get();
  }

  class Iterator
  {
    public:
    explicit Iterator(GttlFastQGenerator* generator, bool end = false)
      : gen(generator), is_end(end)
    {
      if(not end) ++(*this);
    }

    const GttlFastQEntry<buf_size>* operator*() const
    {
      return gen->out;
    }

    const Iterator& operator++()
    {
      if(not gen->advance())
      {
        is_end = true;
      }
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
  GttlLineGenerator<buf_size, true> lg;
};

#endif  // GTTL_FASTQ_GENERATOR_HPP
