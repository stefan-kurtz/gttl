#ifndef GTTL_BINARY_READ_HPP
#define GTTL_BINARY_READ_HPP

#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>

// 4KB is default page size, 8/16 might be more memory friendly,
// but 64KB is the default of std::ifstream. I assume this is reasonable.
// so we set the buffer size to (size_t(1) << 16)

template <typename T,
          size_t buf_size = (size_t(1) << 16)/ sizeof(T)>
class BinaryFileReader
{
  static_assert(buf_size > 0, "Buffer size may not be zero");
  static_assert(std::is_trivially_copyable_v<T>,
                "BinaryFileReader can only work with types that are "
                "trivially copyable.");

  public:
  class Iterator
  {
    private:
    FILE* in_fp;
    T buffer[buf_size];
    size_t buffer_pos;
    size_t buffer_size;

    void fill_buf(void)
    {
      if (in_fp == nullptr)
      {
        buffer_size = 0;
        return;
      }
      buffer_size = std::fread(buffer, sizeof(T), buf_size, in_fp);
      buffer_pos = 0;
      if (buffer_size == 0)
      {
        // In this case we reached EOF
        if (std::fclose(in_fp) != 0)
        {
          throw std::runtime_error("could not close file");
        }
        in_fp = nullptr;
      }
    }

    public:
    explicit Iterator(const std::string &inputfile)
      : in_fp(std::fopen(inputfile.c_str(), "rb"))
      , buffer_pos(0)
      , buffer_size(0)
    {
      if (in_fp == nullptr)
      {
        throw std::runtime_error(std::string("failed to open file: \"") +
                                 inputfile + std::string("\""));
      }
      fill_buf();
    }

    // Objects owning a FILE* should not be copied
    Iterator(const Iterator&) = delete;
    Iterator& operator=(const Iterator&) = delete;

    // Though they may be moved, with the file being closed
    Iterator(Iterator&& other) noexcept
      : in_fp(other.in_fp)
      , buffer_pos(other.buffer_pos)
      , buffer_size(other.buffer_size)
    {
      std::memcpy(buffer, other.buffer, sizeof(buffer));
      other.in_fp = nullptr;
    }

    Iterator& operator=(Iterator&& other) noexcept
    {
      if (this != &other)
      {
        if (in_fp != nullptr)
        {
          if (std::fclose(in_fp) != 0)
          {
            std::cerr << "failed to close file\n";
          }
        }
        in_fp = other.in_fp;
        buffer_pos = other.buffer_pos;
        buffer_size = other.buffer_size;
        std::memcpy(buffer, other.buffer, sizeof(buffer));
        other.in_fp = nullptr;
      }
      return *this;
    }

    ~Iterator(void)
    {
      if (in_fp != nullptr)
      {
        if (std::fclose(in_fp) != 0)
        {
          std::cerr << "failed to close file\n";
        }
      }
    }

    // End Iterator
    Iterator(void)
      : in_fp(nullptr)
      , buffer_pos(0)
      , buffer_size(0)
    { }

    const T& operator*() const
    {
      if (buffer_pos >= buffer_size)
      {
        throw std::out_of_range("Dereferencing end iterator!");
      }
      return buffer[buffer_pos];
    }

    Iterator& operator++()
    {
      if (buffer_pos >= buffer_size)
      {
        throw std::out_of_range("Incrementing end iterator!");
      }
      buffer_pos++;
      if (buffer_pos >= buffer_size)
      {
        fill_buf();
      }
      return *this;
    }

    bool operator==(const Iterator& other) const
    {
      return in_fp == other.in_fp and
             buffer_pos == other.buffer_pos and
             buffer_size == other.buffer_size;
    }

    bool operator!=(const Iterator& other) const
    {
      return not(*this == other);
    }
  };

  private:
  const std::string &inputfile;
  public:
  explicit BinaryFileReader(const std::string &_inputfile)
    : inputfile(_inputfile)
  { }
  [[nodiscard]] Iterator begin(void) const { return Iterator(inputfile); }
  [[nodiscard]] Iterator end(void) const { return Iterator(); }
};
#endif  // GTTL_BINARY_READ_HPP
