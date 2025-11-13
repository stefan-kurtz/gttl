#ifndef GTTL_FASTA_GENERATOR_HPP
#define GTTL_FASTA_GENERATOR_HPP

#include "utilities/gttl_file_open.hpp"
#include "utilities/gttl_line_generator.hpp"
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <ios>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

template <const size_t buf_size = (size_t{1} << size_t{14})>
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

  std::string header;
  std::string sequence;

  GttlFastAEntry()
  {
    header.reserve(buf_size);
    sequence.reserve(buf_size);
  }

  GttlFastAEntry(std::string  _header, std::string  _sequence)
    : header(std::move(_header))
    , sequence(std::move(_sequence))
  {}
};

template <const size_t buf_size = (size_t{1} << size_t{14})>
class GttlFastAGenerator
{
  private:
  // A default output buffer, in case the caller does not provide one
  GttlFastAEntry<buf_size> default_buffer;
  // output-pointer towards the GttlFastAEntry in which the current line should
  // be stored.
  GttlFastAEntry<buf_size>* out;
  // Whether all sequences have been read
  bool is_end;
  // Whether we are still reading the first entry.
  // This flag is required to properly skip the initial FastA-header.
  bool is_first_entry = true;
  // Whether we wish to ignore empty lines
  static constexpr const bool skip_empty_lines = true;
  // The line-generator we use to iterate over the file
  GttlLineGenerator<buf_size, skip_empty_lines> lg;

  public:
  static constexpr const bool is_fastq_generator = false;
  explicit GttlFastAGenerator(const char* file_name,
                              bool _is_end = false)
    : out(&default_buffer)
    , is_end(_is_end)
    , lg(gttl_fp_type_open(file_name, "rb"), &(out->header), _is_end)
  {}

  explicit GttlFastAGenerator(GttlFpType fp,
                              bool _is_end = false)
    : out(&default_buffer)
    , is_end(_is_end)
    , lg(fp, &out->header, _is_end)
  {}

  explicit GttlFastAGenerator(const char* file_name,
                              GttlFastAEntry<buf_size> *_out,
                              bool _is_end = false)
    : out(_out)
    , is_end(_is_end)
    , lg(gttl_fp_type_open(file_name, "rb"), &_out->header, _is_end)
  {}

  explicit GttlFastAGenerator(GttlFpType fp,
                              GttlFastAEntry<buf_size> *_out,
                              bool _is_end = false)
    : out(_out)
    , is_end(_is_end)
    , lg(fp, &_out->header, _is_end)
  {}

  explicit GttlFastAGenerator(const char* _input_string,
                              size_t _string_length,
                              bool _is_end = false)
    : out(&default_buffer)
    , is_end(_is_end)
    , lg(_input_string, _string_length)
  {}

  explicit GttlFastAGenerator(const std::string_view &_input_string,
                              bool _is_end = false)
    : out(&default_buffer)
    , is_end(_is_end)
    , lg(_input_string.data(), _input_string.size())
  {}

  explicit GttlFastAGenerator(const std::vector<std::string>* _file_list,
                              GttlFastAEntry<buf_size> *_out = nullptr,
                              bool _is_end = false)
    : out(_out == nullptr ? &default_buffer : _out)
    , is_end(_is_end)
    , lg(_file_list, &out->header, _is_end)
  {}

  // Delete copy/move constructur & assignment operator
  GttlFastAGenerator(const GttlFastAGenerator&) = delete;
  GttlFastAGenerator& operator=(const GttlFastAGenerator&) = delete;
  GttlFastAGenerator(GttlFastAGenerator&&) = delete;
  GttlFastAGenerator& operator=(GttlFastAGenerator&&) = delete;

  bool advance()
  {
    if(is_end) return false;

    if(out == nullptr)
    {
      lg.set_line_buffer(nullptr);

      if(not std::get<0>(lg.advance())) return false;

      while(true)
      {
        const int ch = lg.getc();
        if(ch == EOF or ch == '>') break;
        lg.advance();
      }

      is_end = (lg.line_number_get() == 0);
      return not is_end;
    }

    if(is_first_entry)
    {
      const int ch = lg.getc();
      if(ch != '>')
      {
        throw std::ios_base::failure(std::string(", line ")
                                     + std::to_string(lg.line_number_get())
                                     + ": corrupted sequence");
      }
      is_first_entry = false;
    }

    out->header.clear();
    out->sequence.clear();
    lg.set_line_buffer(&out->header);

    if(not std::get<0>(lg.advance()))
    {
      is_end = true;
      return false;
    }

    lg.set_line_buffer(&out->sequence);
    while(true)
    {
      const int ch = lg.getc();
      if(ch == EOF or ch == '>') break;

      out->sequence.push_back(static_cast<char>(ch));
      if(not std::get<0>(lg.advance(true))) break;
    }

    if(out->sequence_get().empty() or out->sequence_get()[0] == '>')
    {
      throw std::ios_base::failure(std::string(", line ")
                                   + std::to_string(lg.line_number_get())
                                   + ": corrupted sequence");
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

    const GttlFastAEntry<buf_size>* operator*() const
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
};
#endif // GTTL_FASTA_GENERATOR_HPP
