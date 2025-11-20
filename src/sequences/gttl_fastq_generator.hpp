#ifndef GTTL_FASTQ_GENERATOR_HPP
#define GTTL_FASTQ_GENERATOR_HPP

#include "utilities/gttl_file_open.hpp"
#include "utilities/gttl_line_generator.hpp"
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <ios>
#include <string>
#include <string_view>
#include <vector>

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

  GttlFastQEntry(void)
  {
    header.reserve(buf_size);
    sequence.reserve(buf_size);
    quality.reserve(buf_size);
  }
};

template <const size_t buf_size = (size_t{1} << size_t{14})>
class GttlFastQGenerator
{
  private:
  // A default output buffer, in case the caller does not provide one
  GttlFastQEntry<buf_size> default_buffer;
  // output-pointer towards the GttlFastAEntry in which the current line should
  // be stored.
  GttlFastQEntry<buf_size>* out;
  // Whether all sequences have been read
  bool is_end;
  // Whether we wish to ignore empty lines
  static constexpr const bool skip_empty_lines = true;
  // The line-generator we use to iterate over the file
  GttlLineGenerator<buf_size, skip_empty_lines> lg;

  public:
  static constexpr const bool is_fastq_generator = true;
  explicit GttlFastQGenerator(const char* file_name,
                              bool _is_end = false)
    : out(&default_buffer)
    , is_end(_is_end)
    , lg(gttl_fp_type_open(file_name, "rb"), &(out->header), _is_end)
  { }

  explicit GttlFastQGenerator(GttlFpType fp,
                              bool _is_end = false)
    : out(&default_buffer)
    , is_end(_is_end)
    , lg(fp, &out->header, _is_end)
  { }

  explicit GttlFastQGenerator(const char* file_name,
                              GttlFastQEntry<buf_size> *_out,
                              bool _is_end = false)
    : out(_out)
    , is_end(_is_end)
    , lg(gttl_fp_type_open(file_name, "rb"), _out->header, _is_end)
  { }

  explicit GttlFastQGenerator(GttlFpType fp,
                              GttlFastQEntry<buf_size> *_out,
                              bool _is_end = false)
    : out(_out)
    , is_end(_is_end)
    , lg(fp, _out->header, _is_end)
  { }

  explicit GttlFastQGenerator(const char* _input_string,
                              size_t _string_length,
                              bool _is_end = false)
    : out(&default_buffer)
    , is_end(_is_end)
    , lg(_input_string, _string_length)
  { }

  explicit GttlFastQGenerator(const std::string_view &_input_string,
                              bool _is_end = false)
    : out(&default_buffer)
    , is_end(_is_end)
    , lg(_input_string.data(), _input_string.size())
  { }

  explicit GttlFastQGenerator(const std::vector<std::string>* _file_list,
                              GttlFastQEntry<buf_size> *_out = nullptr,
                              bool _is_end = false)
    : out(_out == nullptr ? &default_buffer : _out)
    , is_end(_is_end)
    , lg(_file_list, &out->header, _is_end)
  { }

  // Delete copy/move constructur & assignment operator
  GttlFastQGenerator(const GttlFastQGenerator&) = delete;
  GttlFastQGenerator& operator = (const GttlFastQGenerator&) = delete;
  GttlFastQGenerator(GttlFastQGenerator&&) = delete;
  GttlFastQGenerator& operator = (GttlFastQGenerator&&) = delete;

  bool advance(void)
  {
    if (is_end)
    {
      return false;
    }
    if (out == nullptr)
    {
      lg.set_line_buffer(nullptr);
      for(size_t line = 0; line < size_t(4); line++)
      {
        lg.advance();
      }
      return true;
    }
    lg.set_line_buffer(&out->header);
    const char cc = lg.getc();
    if (cc == EOF)
    {
      is_end = true;
      return false;
    }
    if (cc != '@')
    {
      throw std::ios_base::failure(std::string(", line ")
                                   + std::to_string(lg.line_number_get())
                                   + ": corrupted sequence");
    }
    lg.advance();
    lg.set_line_buffer(&out->sequence);
    lg.advance();
    lg.set_line_buffer(nullptr);
    lg.advance();
    lg.set_line_buffer(&out->quality);
    return std::get<0>(lg.advance());
  }

  void reset(void)
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
    private:
    GttlFastQGenerator* generator;
    bool is_end;
    public:
    explicit Iterator(GttlFastQGenerator* _generator, bool end = false)
      : generator(_generator)
      , is_end(end)
    {
      if (not end)
      {
        ++(*this);
      }
    }

    const GttlFastQEntry<buf_size>* operator * () const
    {
      return generator->out;
    }

    const Iterator& operator ++ (void)
    {
      assert(not is_end);
      if (not generator->advance())
      {
        is_end = true;
      }
      return *this;
    }

    bool operator == (const Iterator& other) const
    {
      return is_end == other.is_end and generator == other.generator;
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

  Iterator null_iterator(void)
  {
    return Iterator(nullptr, true);
  }
};
#endif  // GTTL_FASTQ_GENERATOR_HPP
