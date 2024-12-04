#ifndef GTTL_LINE_GENERATOR_HPP
#define GTTL_LINE_GENERATOR_HPP

#include <cstdio>
#include <stdexcept>
#include <string>
#include "utilities/gttl_file_open.hpp"


template<size_t buf_size>
class GttlLineGenerator
{
  public:
    class Iterator
    {
      public:
        Iterator(GttlFpType _file, bool _is_end = false)
          : file(_file)
          , is_end(_is_end)
          , buffer(nullptr)
        {
          if(file && !is_end)
          {
            buffer = new char[buf_size];
            advance();
          }
        }

        ~Iterator()
        {
            delete[] buffer;
        }

        Iterator& operator++()
        {
          advance();
          return *this;
        }

        const std::string& operator*() const
        {
          return current_line;
        }

        bool operator==(const Iterator& other) const
        {
          return is_end == other.is_end && file == other.file;
        }

        bool operator!=(const Iterator& other) const
        {
          return !(*this == other);
        }

        private:
        GttlFpType file;
        bool is_end;
        char* buffer;
        std::string current_line;

        void advance()
        {
          if(file && gttl_fp_type_gets(file, buffer, buf_size))
          {
            current_line = std::string(buffer);
            //TODO: Do we want to remove trailing newline characters here?
          }else{
            is_end = true;
          }
        }
    };

    explicit GttlLineGenerator(const char* inputfile)
      : file(gttl_fp_type_open(inputfile, "r"))
    {
      if(!file)
      {
        throw std::runtime_error(std::string("Unable to open file: ")
                                 + inputfile);
      }
    }

    ~GttlLineGenerator()
    {
      if(file)
      {
        gttl_fp_type_close(file);
      }
    }

    // For lvalue-reference iteration
    Iterator begin() &
    {
      return Iterator(file);
    }

    Iterator end() &
    {
      return Iterator(file, true);
    }

    // For rvalue-reference iteration
    Iterator begin() &&
    {
      return Iterator(file);
    }

    Iterator end() &&
    {
      Iterator end_it(file, true);
      file = nullptr;
      own_file = false;
      return end_it;
    }

    void reset(void)
    {
      if(!file)
      {
        throw std::runtime_error("Invalid file handle.");
      }
      gttl_fp_type_reset(file);
    }

    private:
    GttlFpType file;
    bool own_file;
};

#endif // GTTL_LINE_GENERATOR_HPP
