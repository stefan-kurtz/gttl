#ifndef GTTL_BINARY_READ_HPP
#define GTTL_BINARY_READ_HPP

#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>

// 4KB is default page size, 8/16 might be more memory friendly,
// but 64KB is the default of std::ifstream. I assume this is reasonable.
constexpr const size_t default_buffer_bytes = 1U << 16U;

template <typename T,
          const size_t buf_size = (default_buffer_bytes / sizeof(T))>
class BinaryFileIterator
{
 private:
  const char* filename;

 public:
  class Iterator
  {
   public:
    explicit Iterator(const char* filename)
        : file(std::fopen(filename, "rb")), buffer_pos(0), buffer_size(0)
    {
      if (file == nullptr)
      {
        throw std::runtime_error(std::string("Failed to open file: ") +
                                 filename);
      }
      fill_buf();
    }

    // Objects owning a FILE* should not be copied
    Iterator(const Iterator&) = delete;
    Iterator& operator=(const Iterator&) = delete;

    // Though they may be moved
    Iterator(Iterator&& other) noexcept
        : file(other.file),
          buffer_pos(other.buffer_pos),
          buffer_size(other.buffer_size)
    {
      std::memcpy(buffer, other.buffer, sizeof(buffer));
      other.file = nullptr;
    }
    Iterator& operator=(Iterator&& other) noexcept
    {
      if (this != &other)
      {
        if (file != nullptr)
        {
          if (std::fclose(file) != 0)
          {
            std::cerr << "Failed to close file!\n";
          }
        }
        file = other.file;
        buffer_pos = other.buffer_pos;
        buffer_size = other.buffer_size;
        std::memcpy(buffer, other.buffer, sizeof(buffer));
        other.file = nullptr;
      }
      return *this;
    }

    ~Iterator()
    {
      if (file != nullptr)
      {
        if (std::fclose(file) != 0)
        {
          std::cerr << "Failed to close file!\n";
        }
      }
    }

    // End Iterator
    Iterator() : file(nullptr), buffer_pos(0), buffer_size(0) {}

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
      return file == other.file and buffer_pos == other.buffer_pos and
             buffer_size == other.buffer_size;
    }

    bool operator!=(const Iterator& other) const { return !(*this == other); }

   private:
    FILE* file;
    T buffer[buf_size];
    size_t buffer_pos;
    size_t buffer_size;

    void fill_buf()
    {
      if (file == nullptr)
      {
        buffer_size = 0;
        return;
      }
      buffer_size = std::fread(buffer, sizeof(T), buf_size, file);
      buffer_pos = 0;
      if (buffer_size == 0)
      {
        // In this case we reached EOF
        file = nullptr;
      }
    }
  };

  explicit BinaryFileIterator(const char* filename) : filename(filename) {}
  Iterator begin() { return Iterator(filename); }
  Iterator end() { return Iterator(); }
};
#endif  // GTTL_BINARY_READ_HPP
