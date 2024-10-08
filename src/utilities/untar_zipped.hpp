#ifndef UNTAR_ZIPPED_HPP
#define UNTAR_ZIPPED_HPP

#include <cstddef>
#include <cstdlib>
#include <string>
#include <cstdint>
#include "utilities/split_string.hpp"
#include "utilities/has_suffix_or_prefix.hpp"
#include "utilities/popen_reader.hpp"

class DecompressedFile
{
  std::string current_filename;
  size_t current_file_size;
  uint8_t *data_ptr;
  public:
  DecompressedFile(void)
    : current_file_size(0)
    , data_ptr(nullptr)
  {}
  bool is_directory(void) const
  {
    assert(current_filename.size() > 0);
    return current_filename.back() == '/';
  }
  void set(const std::string &_current_filename,
           size_t _current_file_size,
           uint8_t *_data_ptr)
  {
    assert(current_file_size == 0);
    assert(data_ptr == nullptr);
    current_filename = _current_filename;
    current_file_size = _current_file_size;
    data_ptr = _data_ptr;
  }
  const std::string &filename_get(void) const
  {
    return current_filename;
  }
  size_t size(void) const
  {
    return current_file_size;
  }
  uint8_t *data(void) const
  {
    return data_ptr;
  }
  void delete_data(void)
  {
    delete[] data_ptr;
    data_ptr = nullptr;
  }
  void clear(void)
  {
    data_ptr = nullptr;
    current_file_size = 0;
  }
};

class TarReader
{
  PopenReader popen_reader;
  std::string current_filename;
  size_t current_file_size,
         current_file_pos;

  // returns true if next file succesfully setup
  bool next_file(void)
  {
    if (current_file_pos != current_file_size)
    {
      throw std::string("current file ") + current_filename +
            std::string("not read until its end");
    }
    // setup next file
    std::string line;
    while (true)
    {
      const int cc = popen_reader.fgetc_stderr();
      if (cc == '\n')
      {
        break;
      }
      if (cc == EOF)
      {
        return false;
      }
      line.push_back(static_cast<char>(cc));
    }
    std::vector<std::string> line_vector = gttl_split_string(line,' ');
    if (line_vector.size() < 6)
    {
      throw std::string("line \"") + line +
            std::string("\" does not consist of exactly 6 columns");
    }
    if (std::sscanf(line_vector[2].data(),"%zu",&current_file_size) != 1)
    {
      throw std::string("cannot extract byte number from \"") +
            line_vector[2] + std::string("\"");
    }
    line.clear();
    current_filename = std::string(line_vector[5]);
    current_file_pos = 0;
    return true;
  }

  uint8_t *read_from_stdout(size_t count)
  {
    if (count > 0)
    {
      assert(current_file_pos <= current_file_size);
      count = std::min(count, current_file_size - current_file_pos);
      uint8_t *data_ptr = new uint8_t [count];
      const size_t read_bytes = popen_reader.fread_stdout(data_ptr,
                                                          size_t(1), count);
      current_file_pos += read_bytes;
      if (read_bytes != count)
      {
        throw std::string("Error reading file. Couldn't read to finish.");
      }
      return data_ptr;
    }
    return nullptr;
  }
  const char *option_string(const std::string &filename) const
  {
    if (gttl_has_suffix(filename,".tar.bz2"))
    {
      return "-Oxvvjf";
    }
    if (gttl_has_suffix(filename,".tar.gz"))
    {
      return "-Oxvvzf";
    }
    throw std::string("illegal filename ") +
          filename +
          std::string(", can only handle files with suffix "
                      ".tar.gz or .tar.bz2");
  }

  public:
  TarReader(const std::string &filename)
      : popen_reader({"gltar","gtar","tar"},"tar",option_string(filename),
                     filename.c_str())
      , current_file_size(0)
      , current_file_pos(0)
  {
    next_file();
  }
  size_t current_file_size_get(void) const
  {
    return current_file_size;
  }

  class Iterator
  {
    DecompressedFile local_entry;
    TarReader* reader;
    bool exhausted;
    void update_iterator(void)
    {
      local_entry.set(reader->current_filename,
                      reader->current_file_size_get(),
                      reader->read_from_stdout(reader->
                                                 current_file_size_get()));
    }

    public:
    Iterator(TarReader* _reader)
      : reader(_reader)
      , exhausted(false)
    {
      update_iterator();
    }

    Iterator(void)
      : reader(nullptr)
      , exhausted(true)
    { }

    DecompressedFile operator*(void)
    {
      DecompressedFile new_entry{local_entry};
      local_entry.clear();
      return new_entry;
    }

    Iterator& operator++(void)
    {
      exhausted = not reader->next_file();
      if (not exhausted)
      {
        update_iterator();
      }
      return *this;
    }

    bool operator==(const Iterator& other) const
    {
      return exhausted == other.exhausted;
    }

    bool operator!=(const Iterator& other) const
    {
      return not (*this == other);
    }
  };

  Iterator begin(void)
  {
    return Iterator(this);
  }

  Iterator end(void)
  {
    return Iterator();
  }
};
#endif
