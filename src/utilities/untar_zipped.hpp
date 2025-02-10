#ifndef UNTAR_ZIPPED_HPP
#define UNTAR_ZIPPED_HPP

#include <iostream>
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
  const uint8_t *data_ptr;
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
           const uint8_t *_data_ptr)
  {
    assert(current_file_size == 0 and data_ptr == nullptr);
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
  const uint8_t *data(void) const
  {
    return data_ptr;
  }
  void clear(void)
  {
    data_ptr = nullptr;
    current_file_size = 0;
  }
};

class TarReader
{
  uint8_t *data_ptr; /* TarReader keeps its own buffer */
  size_t bytes_allocated,
         current_file_size,
         current_file_pos;
  PopenReader *popen_reader;
  std::string current_filename;
  const bool append_0_byte;

  // returns true if next file succesfully setup
  bool next_file(void)
  {
    if (current_file_pos != current_file_size)
    {
      throw std::string("current file ") + current_filename +
            std::string("not read until its end");
    }
    std::string line;
    while (true)
    {
      const int cc = popen_reader->fgetc_stderr();
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

  const uint8_t *read_from_stdout(size_t count, bool append_0_byte)
  {
    if (count > 0)
    {
      assert(current_file_pos <= current_file_size);
      const size_t bytes2read
        = std::min(count, current_file_size - current_file_pos);
      if (bytes2read > bytes_allocated)
      {
        const size_t extra_byte = append_0_byte ? 1 : 0;
        data_ptr = static_cast<uint8_t *>(realloc(data_ptr,
                                                  sizeof(uint8_t)
                                                  * (bytes2read + extra_byte)));
        bytes_allocated = bytes2read;
      }
      const size_t read_bytes = popen_reader->fread_stdout(data_ptr,
                                                           size_t(1),
                                                           bytes2read);
      if (read_bytes != bytes2read)
      {
        throw std::string("Error reading file: Could not read to finish.");
      }
      if (append_0_byte)
      {
        data_ptr[read_bytes] = '\0';
      }
      current_file_pos += read_bytes;
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
    if (gttl_has_suffix(filename,".tar"))
    {
      return "-Oxvvf";
    }
    throw std::string("illegal filename ") +
          filename +
          std::string(", can only handle files with suffix "
                      ".tar.gz or .tar.bz2 or .tar");
  }

  public:
  TarReader(const std::string &filename,bool with_rapidgzip,bool _append_0_byte)
    : data_ptr(nullptr)
    , bytes_allocated(0)
    , current_file_size(0)
    , current_file_pos(0)
    , popen_reader(nullptr)
    , append_0_byte(_append_0_byte)
  {
    if (with_rapidgzip and not gttl_has_suffix(filename,".tar"))
    {
      popen_reader = new PopenReader({"gtar","tar"},
                                     "tar","-I","rapidgzip","-Oxvvf",
                                     filename.c_str());
    } else
    {
      popen_reader = new PopenReader({"gltar","gtar","tar"},
                                     "tar",option_string(filename),
                                     filename.c_str());
    }
    next_file();
  }
  ~TarReader(void)
  {
    free(data_ptr);
    delete popen_reader;
  }
  size_t current_file_size_get(void) const
  {
    return current_file_size;
  }

  class Iterator
  {
    DecompressedFile local_entry;
    TarReader* reader;
    bool append_0_byte,
         exhausted;
    void update_iterator(void)
    {
      local_entry.set(reader->current_filename,
                      reader->current_file_size_get(),
                      reader->read_from_stdout(reader->
                                                 current_file_size_get(),
                                               append_0_byte));
    }

    public:
    Iterator(TarReader* _reader)
      : reader(_reader)
      , append_0_byte(_reader->append_0_byte)
      , exhausted(false)
    {
      update_iterator();
    }

    Iterator(void)
      : reader(nullptr)
      , append_0_byte(false)
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
